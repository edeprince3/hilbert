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

    // alpha + beta MO transformation matrix
    C_ = (std::shared_ptr<Matrix>)(new Matrix(2*nso_,2*nmo_));
    double ** cp = C_->pointer();
    double ** ca = Ca_->pointer();
    double ** cb = Cb_->pointer();

    for (int mu = 0; mu < nso_; mu++) {
        long int count = 0;
        for (long int i = 0; i < nalpha_; i++) {
            cp[mu][count++] = ca[mu][i];
        }
        for (long int i = 0; i < nbeta_; i++) {
            cp[mu+nso_][count++] = cb[mu][i];
        }
        for (long int i = nalpha_; i < nmo_; i++) {
            cp[mu][count++] = ca[mu][i];
        }
        for (long int i = nbeta_; i < nmo_; i++) {
            cp[mu+nso_][count++] = cb[mu][i];
        }
    }

    std::shared_ptr<Matrix> F (new Matrix(2*nso_,2*nso_));
    double ** fp = F->pointer();
    double ** fa = Fa_->pointer();
    double ** fb = Fb_->pointer();
    for (int mu = 0; mu < nso_; mu++) {
        for (int nu = 0; nu < nso_; nu++) {
            fp[mu][nu] = fa[mu][nu];
            fp[mu+nso_][nu+nso_] = fb[mu][nu];
        }
    }

    // get MO-basis quantities:
    F->transform(C_);

    // orbital energies
    epsilon_ = (double*)malloc(2*nso_*sizeof(double));
    memset((void*)epsilon_,'\0',2*nso_*sizeof(double));

    double ** eps = F->pointer();
    for (long int i = 0; i < 2*nso_; i++) {
        epsilon_[i] = eps[i][i];
    }

    build_mo_eris();

    // allocate memory for amplitudes, residual, and temporary buffer

    long int o = nalpha_ + nbeta_;
    long int v = (nmo_-nalpha_) + (nmo_-nbeta_);

    tamps_ = (double*)malloc((o*v + o*o*v*v)*sizeof(double));
    t2_    = tamps_;
    t1_    = tamps_ + o*o*v*v;

    residual_ = (double*)malloc((o*v + o*o*v*v)*sizeof(double));
    r2_       = residual_;
    r1_       = residual_ + o*o*v*v;

    memset((void*)tamps_,'\0',(o*o*v*v+o*v)*sizeof(double));
    memset((void*)residual_,'\0',(o*o*v*v+o*v)*sizeof(double));

    // temporary storage ... reduce later

    tmp1_ = (double*)malloc(o*v*v*v*sizeof(double));
    tmp2_ = (double*)malloc(o*o*v*v*sizeof(double));
    tmp3_ = (double*)malloc(o*o*v*v*sizeof(double));
    memset((void*)tmp1_,'\0',o*v*v*v*sizeof(double));
    memset((void*)tmp2_,'\0',o*o*v*v*sizeof(double));
    memset((void*)tmp3_,'\0',o*o*v*v*sizeof(double));

    // initialize diis solver
    diis = (std::shared_ptr<DIIS>)(new DIIS(o*o*v*v + o*v));

}

// TODO: this is so wasteful. need to
// 1. fold t1 into 3-index integrals
// 2. don't build full four-index tensor
// 3. don't build (ac|bd) 

void PolaritonicUCCSD::build_mo_eris() {

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    // TODO: use DF_BASIS_CC

    // get auxiliary basis:
    std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

    // total number of auxiliary basis functions
    nQ_ = auxiliary->nbf();

    // three-index integrals
    std::shared_ptr<DFTensor> DFa (new DFTensor(primary,auxiliary,Ca_,nalpha_,nso_-nalpha_,nalpha_,nso_-nalpha_,options_));

    std::shared_ptr<Matrix> tmp_so = DFa->Qso();
    double ** tmp_so_p = tmp_so->pointer();

    std::shared_ptr<Matrix> Qso (new Matrix(nQ_,2*nso_*2*nso_));
    double ** qso_p = Qso->pointer();
    for (int Q = 0; Q < nQ_; Q++) {
        for (int mu = 0; mu < nso_; mu++) {
            for (int nu = 0; nu < nso_; nu++) {
                qso_p[Q][(mu     )*2*nso_+(nu     )] = tmp_so_p[Q][mu*nso_+nu];
                qso_p[Q][(mu+nso_)*2*nso_+(nu+nso_)] = tmp_so_p[Q][mu*nso_+nu];
            }
        }
    }

    // AO->MO transformation

    long int n  = 2L*(long int)nmo_;
    long int ns = 2L*(long int)nso_;

    double * Qmo = (double*)malloc(nQ_*n*n*sizeof(double));
    memset((void*)Qmo,'\0',nQ_*n*n*sizeof(double));

    double * tmp = (double*)malloc(nQ_*n*n*sizeof(double));
    memset((void*)tmp,'\0',nQ_*n*n*sizeof(double));

    // I(Q,mu,p) = C(nu,p) Qso(Q,mu,nu)
    F_DGEMM('n','n',n,ns*nQ_,ns,1.0,&(C_->pointer()[0][0]),n,&(qso_p[0][0]),ns,0.0,Qmo,n);
    for (int Q = 0; Q < nQ_; Q++) {
        for (int p = 0; p < n; p++) {
            for (int mu = 0; mu < ns; mu++) {
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

    // unpack different classes of eris

    long int o = nalpha_ + nbeta_;
    long int v = (nmo_-nalpha_) + (nmo_-nbeta_);

    // <ij||kl>
    eri_ijkl_ = (double*)malloc(o*o*o*o*sizeof(double));
    memset((void*)eri_ijkl_,'\0',o*o*o*o*sizeof(double));

    for (long int i = 0; i < o; i++) {
        for (long int j = 0; j < o; j++) {
            for (long int k = 0; k < o; k++) {
                for (long int l = 0; l < o; l++) {
                    long int ikjl = i*n*n*n+k*n*n+j*n+l;
                    long int iljk = i*n*n*n+l*n*n+j*n+k;
                    eri_ijkl_[i*o*o*o+j*o*o+k*o+l] = eri[ikjl] - eri[iljk];
                }
            }
        }
    }

    // <ab||cd>
    eri_abcd_ = (double*)malloc(v*v*v*v*sizeof(double));
    memset((void*)eri_abcd_,'\0',v*v*v*v*sizeof(double));

    for (long int a = 0; a < v; a++) {
        for (long int b = 0; b < v; b++) {
            for (long int c = 0; c < v; c++) {
                for (long int d = 0; d < v; d++) {
                    long int acbd = (a+o)*n*n*n+(c+o)*n*n+(b+o)*n+(d+o);
                    long int adbc = (a+o)*n*n*n+(d+o)*n*n+(b+o)*n+(c+o);
                    eri_abcd_[a*v*v*v+b*v*v+c*v+d] = eri[acbd] - eri[adbc];
                }
            }
        }
    }

    // <ij||ab>
    eri_ijab_ = (double*)malloc(o*o*v*v*sizeof(double));
    memset((void*)eri_ijab_,'\0',o*o*v*v*sizeof(double));

    for (long int i = 0; i < o; i++) {
        for (long int j = 0; j < o; j++) {
            for (long int a = 0; a < v; a++) {
                for (long int b = 0; b < v; b++) {
                    long int iajb = i*n*n*n+(a+o)*n*n+j*n+(b+o);
                    long int ibja = i*n*n*n+(b+o)*n*n+j*n+(a+o);
                    eri_ijab_[i*o*v*v+j*v*v+a*v+b] = eri[iajb] - eri[ibja];
                }
            }
        }
    }

    // <ia||jb>
    eri_iajb_ = (double*)malloc(o*o*v*v*sizeof(double));
    memset((void*)eri_iajb_,'\0',o*o*v*v*sizeof(double));
    for (long int a = 0; a < v; a++) {
        for (long int b = 0; b < v; b++) {
            for (long int i = 0; i < o; i++) {
                for (long int j = 0; j < o; j++) {
                    long int ijab = i*n*n*n+j*n*n+(a+o)*n+(b+o);
                    long int ibaj = i*n*n*n+(b+o)*n*n+(a+o)*n+j;
                    eri_iajb_[i*o*v*v+a*o*v+j*v+b] = eri[ijab] - eri[ibaj];
                }
            }
        }
    }

    // <ia||jk>
    eri_iajk_ = (double*)malloc(o*o*o*v*sizeof(double));
    memset((void*)eri_iajk_,'\0',o*o*o*v*sizeof(double));
    for (long int i = 0; i < o; i++) {
        for (long int a = 0; a < v; a++) {
            for (long int j = 0; j < o; j++) {
                for (long int k = 0; k < o; k++) {
                    long int ijak = i*n*n*n+j*n*n+(a+o)*n+k;
                    long int ikaj = i*n*n*n+k*n*n+(a+o)*n+j;
                    eri_iajk_[i*o*o*v+a*o*o+j*o+k] = eri[ijak] - eri[ikaj];
                }
            }
        }
    }

    // <ab||ci>
    eri_aibc_ = (double*)malloc(o*v*v*v*sizeof(double));
    memset((void*)eri_aibc_,'\0',o*v*v*v*sizeof(double));

    for (long int a = 0; a < v; a++) {
        for (long int b = 0; b < v; b++) {
            for (long int c = 0; c < v; c++) {
                for (long int i = 0; i < o; i++) {
                    long int abic = (a+o)*n*n*n+(b+o)*n*n+i*n+(c+o);
                    long int acib = (a+o)*n*n*n+(c+o)*n*n+i*n+(b+o);
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
    long int maxiter          = options_.get_int("MAXITER");

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

    long int o = nalpha_ + nbeta_;
    long int v = 2 * nalpha_ - o;

    double ec = 0.0;

    long int iter = 0;
    do {

        e_last = energy_ + ec;

        // build residual 
        residual();

        // update amplitudes 
        tnorm = update_amplitudes();

        // evaluate correlation energy
        ec = correlation_energy();

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

void PolaritonicUCCSD::residual() {

    long int o = nalpha_ + nbeta_;
    long int v = 2 * nmo_ - o;

    memset((void*)residual_,'\0',(o*o*v*v+o*v)*sizeof(double));

    // singles
    // 
    // pdaggerq's unfactorized output:
    // 
    // < 0 | m* e e(-T) H e(T) | 0> :
    // 
    //     + 1.00000 F(e,m) 
    //     - 1.00000 F(i,m) t1(e,i) 
    //     + 1.00000 F(e,a) t1(a,m) 
    //     + 1.00000 F(i,a) t2(a,e,i,m) 
    //     - 1.00000 F(i,a) t1(a,m) t1(e,i) 
    //     - 1.00000 <i,e||m,a> t1(a,i) 
    //     + 0.50000 <i,j||m,a> ( t2(a,e,i,j) + 2 t1(a,i) t1(e,j) )
    //     - 0.50000 <e,i||a,b> ( t2(a,b,i,m) + 2 t1(a,i) t1(b,m) )
    //     - 1.00000 <i,j||a,b> t1(a,i) ( t2(b,e,m,j) + t1(b,m) t1(e,j) )
    //     + 0.50000 <i,j||a,b> t1(a,m) t2(b,e,i,j) 
    //     + 0.50000 <i,j||a,b> t1(e,i) t2(a,b,j,m) 

    // optimized:

    // - 1.00000 <i,e||m,a> t1(a,i) 

#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int m = 0; m < o; m++) {
            double dum = 0.0;
            for (long int a = 0; a < v; a++) {
                for (long int i = 0; i < o; i++) {
                    dum -= eri_iajb_[i*o*v*v+e*o*v+m*v+a] * t1_[a*o+i];
                }
            }
            r1_[e*o+m] = dum;
        }
    }

    // + 0.50000 <i,j||m,a> ( t2(a,e,i,j) + 2 t1(a,i) t1(e,j) )

    // v'(m,i,j,a) = <i,j||m,a>
#pragma omp parallel for schedule(static)
    for (long int m = 0; m < o; m++) {
        for (long int i = 0; i < o; i++) {
            for (long int j = 0; j < o; j++) {
                for (long int a = 0; a < v; a++) {
                    tmp1_[m*o*o*v+i*o*v+j*v+a] = eri_iajk_[m*o*o*v+a*o*o+i*o+j];
                }
            }
        }
    }

    // t'(i,j,a,e) = t2(a,e,i,j) + 2 t1(a,i) t1(e,j)
#pragma omp parallel for schedule(static)
    for (long int i = 0; i < o; i++) {
        for (long int j = 0; j < o; j++) {
            for (long int a = 0; a < v; a++) {
                for (long int e = 0; e < v; e++) {
                    tmp2_[i*o*v*v+j*v*v+a*v+e] = t2_[a*o*o*v+e*o*o+i*o+j] + 2.0 * t1_[a*o+i] * t1_[e*o+j];
                }
            }
        }
    }

    // r(e,m) = v'(m,i,j,a) t'(i,j,a,e)
    F_DGEMM('t','t',o,v,o*o*v,0.5,tmp1_,o*o*v,tmp2_,v,1.0,r1_,o);

    // - 0.50000 <e,i||a,b> ( t2(a,b,i,m) + 2 t1(a,i) t1(b,m) )

    // t'(i,a,b,m) = t2(a,b,i,m) + 2 t1(a,i) t1(b,m)
#pragma omp parallel for schedule(static)
    for (long int i = 0; i < o; i++) {
        for (long int a = 0; a < v; a++) {
            for (long int b = 0; b < v; b++) {
                for (long int m = 0; m < o; m++) {
                    tmp1_[i*o*v*v+a*o*v+b*o+m] = t2_[a*o*o*v+b*o*o+i*o+m] + 2.0 * t1_[a*o+i] * t1_[b*o+m];
                }   
            }   
        }   
    }

    // r(e,m) = t'(i,a,b,m) <e,i||a,b>
    F_DGEMM('n','n',o,v,o*v*v,-0.5,tmp1_,o,eri_aibc_,o*v*v,1.0,r1_,o);

    // - 1.00000 <i,j||a,b> t1(a,i) ( t2(b,e,m,j) + t1(b,m) t1(e,j) )
            
    // I(j,b) = <i,j||a,b> t1(a,i)
#pragma omp parallel for schedule(static)
    for (long int j = 0; j < o; j++) {
        for (long int b = 0; b < v; b++) {
            double dum = 0.0;
            for (long int a = 0; a < v; a++) {
                for (long int i = 0; i < o; i++) {
                    dum += eri_ijab_[i*o*v*v+j*v*v+a*v+b] * t1_[a*o+i];
                }
            }   
            tmp1_[j*v+b] = dum;
        }       
    }       

    // r(e,m) = - I(j,b) ( t2(b,e,m,j) + t1(b,m) t1(e,j) )
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int m = 0; m < o; m++) {
            double dum = 0.0;
            for (long int j = 0; j < o; j++) {
                for (long int b = 0; b < v; b++) {
                    dum += tmp1_[j*v+b] * ( t2_[b*o*o*v+e*o*o+m*o+j] + t1_[b*o+m] * t1_[e*o+j] );
                }
            }
            r1_[e*o+m] -= dum;
        }
    }

    // + 0.50000 <i,j||a,b> t1(a,m) t2(b,e,i,j) 

    // v'(i,j,b,a) = <i,j||a,b>
#pragma omp parallel for schedule(static)
    for (long int i = 0; i < o; i++) {
        for (long int j = 0; j < o; j++) {
            for (long int b = 0; b < v; b++) {
                for (long int a = 0; a < v; a++) {
                    tmp1_[i*o*v*v+j*v*v+b*v+a] = eri_ijab_[i*o*v*v+j*v*v+a*v+b];
                }
            }
        }
    }

    // I(i,j,b,m) = t(a,m) v'(i,j,b,a)
    F_DGEMM('n','n',o,o*o*v,v,1.0,t1_,o,tmp1_,v,0.0,tmp2_,o);

    // t'(e,i,j,b) = t2(b,e,i,j)
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) { 
        for (long int i = 0; i < o; i++) {
            for (long int j = 0; j < o; j++) {
                for (long int b = 0; b < v; b++) {
                    tmp1_[e*o*o*v+i*o*v+j*v+b] = t2_[b*o*o*v+e*o*o+i*o+j];
                }
            }
        }
    }

    // r(e,m) = 0.5 I(i,j,b,m) t'(e,i,j,b)
    F_DGEMM('n','n',o,v,o*o*v,0.5,tmp2_,o,tmp1_,o*o*v,1.0,r1_,o);

    // + 0.50000 <i,j||a,b> t1(e,i) t2(a,b,j,m) 

    // t'(j,a,b,m) = t2(a,b,j,m) 
#pragma omp parallel for schedule(static)
    for (long int j = 0; j < o; j++) { 
        for (long int a = 0; a < v; a++) {
            for (long int b = 0; b < v; b++) {
                for (long int m = 0; m < o; m++) {
                    tmp1_[j*o*v*v+a*o*v+b*o+m] = t2_[a*o*o*v+b*o*o+j*o+m];
                }
            }
        }
    }

    // I(i,m) = t'(j,a,b,m) <i,j||a,b>
    F_DGEMM('n','n',o,o,o*v*v,1.0,tmp1_,o,eri_ijab_,o*v*v,0.0,tmp2_,o);

    // r(e,m) = 0.5 I(i,m) t1(e,i)
    F_DGEMM('n','n',o,v,o,0.5,tmp2_,o,t1_,o,1.0,r1_,o);

    // doubles:
    //
    // pdaggerq's unoptimized output:
    //
    //< 0 | m* n* f e e(-T) H e(T) | 0> :
    //
    //
    // fully-contracted strings:
    //     + 1.00000 h(i,n) t2(e,f,i,m) 
    //     - 1.00000 h(i,m) t2(e,f,i,n) 
    //     + 1.00000 h(e,a) t2(a,f,m,n) 
    //     - 1.00000 h(f,a) t2(a,e,m,n) 
    //     - 1.00000 h(i,a) t1(a,n) t2(e,f,m,i) 
    //     + 1.00000 h(i,a) t1(a,m) t2(e,f,n,i) 
    //     - 1.00000 h(i,a) t1(e,i) t2(a,f,m,n) 
    //     + 1.00000 h(i,a) t1(f,i) t2(a,e,m,n) 
    //     + 1.00000 <e,f||m,n> 
    //     + 1.00000 <i,e||m,n> t1(f,i) 
    //     - 1.00000 <i,f||m,n> t1(e,i) 
    //     + 1.00000 <e,f||a,n> t1(a,m) 
    //     - 1.00000 <e,f||a,m> t1(a,n) 
    //     + 1.00000 <i,j||i,n> t2(e,f,j,m) 
    //     - 1.00000 <i,j||i,m> t2(e,f,j,n) 
    //     + 0.50000 <i,j||m,n> t2(e,f,i,j) 
    //     + 1.00000 <i,e||i,a> t2(a,f,m,n) 
    //     - 1.00000 <i,e||a,n> t2(a,f,i,m) 
    //     + 1.00000 <i,e||a,m> t2(a,f,i,n) 
    //     - 1.00000 <i,f||i,a> t2(a,e,m,n) 
    //     + 1.00000 <i,f||a,n> t2(a,e,i,m) 
    //     - 1.00000 <i,f||a,m> t2(a,e,i,n) 
    //     + 0.50000 <e,f||a,b> t2(a,b,m,n) 
    //     + 1.00000 <i,j||m,n> t1(e,i) t1(f,j) 
    //     + 1.00000 <i,e||a,n> t1(a,m) t1(f,i) 
    //     - 1.00000 <i,e||a,m> t1(a,n) t1(f,i) 
    //     - 1.00000 <i,f||a,n> t1(a,m) t1(e,i) 
    //     + 1.00000 <i,f||a,m> t1(a,n) t1(e,i) 
    //     - 1.00000 <e,f||a,b> t1(a,n) t1(b,m) 
    //     - 1.00000 <i,j||i,a> t1(a,n) t2(e,f,m,j) 
    //     + 1.00000 <i,j||i,a> t1(a,m) t2(e,f,n,j) 
    //     - 1.00000 <i,j||a,j> t1(e,i) t2(a,f,m,n) 
    //     + 1.00000 <i,j||a,j> t1(f,i) t2(a,e,m,n) 
    //     + 1.00000 <i,j||a,n> t1(a,i) t2(e,f,j,m) 
    //     + 0.50000 <i,j||a,n> t1(a,m) t2(e,f,i,j) 
    //     - 1.00000 <i,j||a,n> t1(e,i) t2(a,f,j,m) 
    //     + 1.00000 <i,j||a,n> t1(f,i) t2(a,e,j,m) 
    //     - 1.00000 <i,j||a,m> t1(a,i) t2(e,f,j,n) 
    //     - 0.50000 <i,j||a,m> t1(a,n) t2(e,f,i,j) 
    //     + 1.00000 <i,j||a,m> t1(e,i) t2(a,f,j,n) 
    //     - 1.00000 <i,j||a,m> t1(f,i) t2(a,e,j,n) 
    //     + 1.00000 <i,e||a,b> t1(a,i) t2(b,f,m,n) 
    //     - 1.00000 <i,e||a,b> t1(a,n) t2(b,f,m,i) 
    //     + 1.00000 <i,e||a,b> t1(a,m) t2(b,f,n,i) 
    //     + 0.50000 <i,e||a,b> t1(f,i) t2(a,b,m,n) 
    //     - 1.00000 <i,f||a,b> t1(a,i) t2(b,e,m,n) 
    //     + 1.00000 <i,f||a,b> t1(a,n) t2(b,e,m,i) 
    //     - 1.00000 <i,f||a,b> t1(a,m) t2(b,e,n,i) 
    //     - 0.50000 <i,f||a,b> t1(e,i) t2(a,b,m,n) 
    //     - 0.50000 <i,j||a,b> t2(a,b,n,j) t2(e,f,m,i) 
    //     + 0.50000 <i,j||a,b> t2(a,b,m,j) t2(e,f,n,i) 
    //     + 0.25000 <i,j||a,b> t2(a,b,m,n) t2(e,f,i,j) 
    //     - 0.50000 <i,j||a,b> t2(a,e,i,j) t2(b,f,m,n) 
    //     + 1.00000 <i,j||a,b> t2(a,e,n,j) t2(b,f,m,i) 
    //     - 1.00000 <i,j||a,b> t2(a,e,m,j) t2(b,f,n,i) 
    //     - 0.50000 <i,j||a,b> t2(a,e,m,n) t2(b,f,i,j) 
    //     - 1.00000 <i,j||a,n> t1(a,m) t1(e,j) t1(f,i) 
    //     + 1.00000 <i,j||a,m> t1(a,n) t1(e,j) t1(f,i) 
    //     - 1.00000 <i,e||a,b> t1(a,n) t1(b,m) t1(f,i) 
    //     + 1.00000 <i,f||a,b> t1(a,n) t1(b,m) t1(e,i) 
    //     - 1.00000 <i,j||a,b> t1(a,i) t1(b,n) t2(e,f,m,j) 
    //     + 1.00000 <i,j||a,b> t1(a,i) t1(b,m) t2(e,f,n,j) 
    //     - 1.00000 <i,j||a,b> t1(a,i) t1(e,j) t2(b,f,m,n) 
    //     + 1.00000 <i,j||a,b> t1(a,i) t1(f,j) t2(b,e,m,n) 
    //     - 0.50000 <i,j||a,b> t1(a,n) t1(b,m) t2(e,f,i,j) 
    //     + 1.00000 <i,j||a,b> t1(a,n) t1(e,j) t2(b,f,m,i) 
    //     - 1.00000 <i,j||a,b> t1(a,n) t1(f,j) t2(b,e,m,i) 
    //     - 1.00000 <i,j||a,b> t1(a,m) t1(e,j) t2(b,f,n,i) 
    //     + 1.00000 <i,j||a,b> t1(a,m) t1(f,j) t2(b,e,n,i) 
    //     + 0.50000 <i,j||a,b> t1(e,i) t1(f,j) t2(a,b,m,n) 
    //     - 1.00000 <i,j||a,b> t1(a,n) t1(b,m) t1(e,i) t1(f,j) 

    // optimized

    // + 1.00000 <e,f||m,n> 
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[e*o*o*v+f*o*o+m*o+n] = eri_ijab_[m*o*v*v+n*v*v+e*v+f];
                }
            }
        }
    }

    // + 1.00000 <i,e||m,n> t1(f,i) 
    // - 1.00000 <i,f||m,n> t1(e,i) 
    F_DGEMM('n','n',o*o*v,v,o,1.0,eri_iajk_,o*o*v,t1_,o,0.0,tmp1_,o*o*v);
    C_DAXPY(o*o*v*v,-1.0,tmp1_,1,r2_,1);
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[e*o*o*v+f*o*o+m*o+n] += tmp1_[f*o*o*v+e*o*o+m*o+n];
                }
            }
        }
    }


    // - 0.50000 <i,j||n,a> t1(a,m) t2(e,f,i,j) 
    // + 0.50000 <i,j||m,a> t1(a,n) t2(e,f,i,j) 

    // v'(a,n,i,j) = <n,a||i,j>
#pragma omp parallel for schedule(static)
    for (long int a = 0; a < v; a++) {
        for (long int n = 0; n < o; n++) {
            for (long int i = 0; i < o; i++) {
                for (long int j = 0; j < o; j++) {
                    tmp1_[a*o*o*o+n*o*o+i*o+j] = eri_iajk_[n*o*o*v+a*o*o+i*o+j];
                }
            }
        }
    }

    // I(m,n,i,j) = v'(a,n,i,j) t1(a,m)
    F_DGEMM('n','t',o*o*o,o,v,1.0,tmp1_,o*o*o,t1_,o,0.0,tmp2_,o*o*o);

    // r'(e,f,m,n) = 0.5 I(m,n,i,j) t2(e,f,i,j)
    F_DGEMM('t','n',o*o,v*v,o*o,0.5,tmp2_,o*o,t2_,o*o,0.0,tmp1_,o*o);

    C_DAXPY(o*o*v*v,-1.0,tmp1_,1,r2_,1);
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[e*o*o*v+f*o*o+m*o+n] += tmp1_[e*o*o*v+f*o*o+n*o+m];
                }
            }
        }
    }

    // + 1.00000 <i,j||n,a> t1(a,i) t2(e,f,m,j) 
    // - 1.00000 <i,j||m,a> t1(a,i) t2(e,f,n,j)

    // I(j,n) = <i,j||n,a> t1(a,i)
#pragma omp parallel for schedule(static)
    for (long int j = 0; j < o; j++) {
        for (long int n = 0; n < o; n++) {
            double dum = 0.0;
            for (long int a = 0; a < v; a++) {
                for (long int i = 0; i < o; i++) {
                    dum += eri_iajk_[n*o*o*v+a*o*o+i*o+j] * t1_[a*o+i];
                }
            }
            tmp1_[j*o+n] = dum;
        }
    }

    // r'(e,f,m,n) = I(j,n) t2(e,f,m,j)
    F_DGEMM('n','n',o,o*v*v,o,1.0,tmp1_,o,t2_,o,0.0,tmp2_,o);
    C_DAXPY(o*o*v*v,1.0,tmp2_,1,r2_,1);
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[e*o*o*v+f*o*o+m*o+n] -= tmp2_[e*o*o*v+f*o*o+n*o+m];
                }
            }
        }
    }

    // - 1.00000 <i,j||n,a> t1(f,i) t2(a,e,j,m)
    // + 1.00000 <i,j||m,a> t1(f,i) t2(a,e,j,n)
    // + 1.00000 <i,j||n,a> t1(e,i) t2(a,f,j,m) 
    // - 1.00000 <i,j||m,a> t1(e,i) t2(a,f,j,n) 

    // v'(n,i,a,j) = <i,j||n,a>
#pragma omp parallel for schedule(static)
    for (long int n = 0; n < o; n++) {
        for (long int i = 0; i < o; i++) {
            for (long int a = 0; a < v; a++) {
                for (long int j = 0; j < o; j++) {
                    tmp1_[n*o*o*v+i*o*v+a*o+j] = eri_iajk_[n*o*o*v+a*o*o+i*o+j];
                }
            }
        }
    }

    // t'(e,m,a,j) = t2(a,e,j,m)
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int m = 0; m < o; m++) {
            for (long int a = 0; a < v; a++) {
                for (long int j = 0; j < o; j++) {
                    tmp2_[e*o*o*v+m*o*v+a*o+j] = t2_[a*o*o*v+e*o*o+j*o+m];
                }
            }
        }
    }

    // I(e,m,n,i) = v'(n,i,a,j) t'(e,m,a,j)
    F_DGEMM('t','n',o*o,o*v,o*v,1.0,tmp1_,o*v,tmp2_,o*v,0.0,tmp3_,o*o);

    // r'(f,e,m,n) = I(e,m,n,i) t1(f,i)
    F_DGEMM('t','n',o*o*v,v,o,1.0,tmp3_,o,t1_,o,0.0,tmp1_,o*o*v);

#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[e*o*o*v+f*o*o+m*o+n] -= tmp1_[f*o*o*v+e*o*o+m*o+n];
                    r2_[e*o*o*v+f*o*o+m*o+n] += tmp1_[f*o*o*v+e*o*o+n*o+m];
                    r2_[e*o*o*v+f*o*o+m*o+n] += tmp1_[e*o*o*v+f*o*o+m*o+n];
                    r2_[e*o*o*v+f*o*o+m*o+n] -= tmp1_[e*o*o*v+f*o*o+n*o+m];
                }
            }
        }
    }

    //     + 1.00000 <i,j||n,a> t1(f,i) t1(a,m) t1(e,j)
    //     - 1.00000 <i,j||m,a> t1(f,i) t1(a,n) t1(e,j)

    // v'(a,n,i,j) = <na||ij>
#pragma omp parallel for schedule(static)
    for (long int a = 0; a < v; a++) {
        for (long int n = 0; n < o; n++) {
            for (long int i = 0; i < o; i++) {
                for (long int j = 0; j < o; j++) {
                    tmp1_[a*o*o*o+n*o*o+i*o+j] = eri_iajk_[n*o*o*v+a*o*o+i*o+j];
                }
            }
        }
    }

    // I(m,n,i,j) = v'(a,n,i,j) t1(a,m)
    F_DGEMM('n','t',o*o*o,o,v,1.0,tmp1_,o*o*o,t1_,o,0.0,tmp2_,o*o*o);

    // I'(e,m,n,i) = I(m,n,i,j) t1(e,j)
    F_DGEMM('t','n',o*o*o,v,o,1.0,tmp2_,o,t1_,o,0.0,tmp1_,o*o*o);

    // r'(f,e,m,n) = I'(e,m,n,i) t1(f,i)
    F_DGEMM('t','n',o*o*v,v,o,1.0,tmp1_,o,t1_,o,0.0,tmp2_,o*o*v);

#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[e*o*o*v+f*o*o+m*o+n] += tmp2_[f*o*o*v+e*o*o+m*o+n];
                    r2_[e*o*o*v+f*o*o+m*o+n] -= tmp2_[f*o*o*v+e*o*o+n*o+m];
                }
            }
        }
    }

    // + 1.00000 <e,f||a,n> t1(a,m) 
    // - 1.00000 <e,f||a,m> t1(a,n) 
    F_DGEMM('n','t',o*v*v,o,v,1.0,eri_aibc_,o*v*v,t1_,o,0.0,tmp1_,o*v*v);
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[e*o*o*v+f*o*o+m*o+n] += tmp1_[m*o*v*v+n*v*v+e*v+f];
                    r2_[e*o*o*v+f*o*o+m*o+n] -= tmp1_[n*o*v*v+m*v*v+e*v+f];
                }
            }
        }
    }

    // t'(a,b,i,j) = t2(a,b,i,j) + 2 t1(a,i) t1(b,j)
    C_DCOPY(o*o*v*v,t2_,1,tmp1_,1);
#pragma omp parallel for schedule(static)
    for (long int a = 0; a < v; a++) {
        for (long int b = 0; b < v; b++) {
            for (long int i = 0; i < o; i++) {
                for (long int j = 0; j < o; j++) {
                    tmp1_[a*o*o*v+b*o*o+i*o+j] += 2.0 * t1_[a*o+i] * t1_[b*o+j];
                }
            }
        }
    }

    // + 0.50000 <i,j||m,n> t2(e,f,i,j) + 2 t1(e,i) t1(f,j) 
    F_DGEMM('n','n',o*o,v*v,o*o,0.5,eri_ijkl_,o*o,tmp1_,o*o,1.0,residual_,o*o);

    // + 0.50000 <e,f||a,b> t2(a,b,m,n) + 2 t1(a,m) t1(b,n) 
    F_DGEMM('n','n',o*o,v*v,v*v,0.5,tmp1_,o*o,eri_abcd_,v*v,1.0,residual_,o*o);

    // - 1.00000 P(e,f) P(m,n) <i,e||n,a> t2(a,f,m,i) + t1(a,m) t1(f,i) 

    // t'(a,i,f,m) = t2(a,f,m,i) + t1(a,m) t1(f,i) 
#pragma omp parallel for schedule(static)
    for (long int a = 0; a < v; a++) {
        for (long int i = 0; i < o; i++) {
            for (long int f = 0; f < v; f++) {
                for (long int m = 0; m < o; m++) {
                    tmp1_[a*o*o*v+i*o*v+f*o+m] = t2_[a*o*o*v+f*o*o+m*o+i] + t1_[a*o+m] * t1_[f*o+i];
                }
            }
        }
    }

    // v'(a,i,e,n) = <ie||na>
#pragma omp parallel for schedule(static)
    for (long int a = 0; a < v; a++) {
        for (long int i = 0; i < o; i++) {
            for (long int e = 0; e < v; e++) {
                for (long int n = 0; n < o; n++) {
                    tmp2_[a*o*o*v+i*o*v+e*o+n] = eri_iajb_[i*o*v*v+e*o*v+n*v+a];
                }
            }
        }
    }

    // I(e,n,f,m) = t'(a,i,f,m) v'(a,i,e,n) 
    F_DGEMM('n','t',o*v,o*v,o*v,1.0,tmp1_,o*v,tmp2_,o*v,0.0,tmp3_,o*v);
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {

                    double dum = 0.0;
                    dum -= tmp3_[e*o*o*v+n*o*v+f*o+m];
                    dum += tmp3_[e*o*o*v+m*o*v+f*o+n];
                    dum += tmp3_[f*o*o*v+n*o*v+e*o+m];
                    dum -= tmp3_[f*o*o*v+m*o*v+e*o+n];

                    r2_[e*o*o*v+f*o*o+m*o+n] += dum;
                }
            }
        }
    }

    // - 0.50000 <e,i||a,b> t1(f,i) t2(a,b,m,n) 
    // + 0.50000 <f,i||a,b> t1(e,i) t2(a,b,m,n) 

    // I(m,n,e,i) = <e,i||a,b> t2(a,b,m,n)
    F_DGEMM('t','t',o*v,o*o,v*v,0.5,eri_aibc_,v*v,t2_,o*o,0.0,tmp1_,o*v);

    // R(m,n,e,f) = t1(f,i) I(m,n,e,i)
    F_DGEMM('t','n',v,o*o*v,o,1.0,t1_,o,tmp1_,o,0.0,tmp2_,v);
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[e*o*o*v+f*o*o+m*o+n] -= tmp2_[m*o*v*v+n*v*v+e*v+f];
                    r2_[e*o*o*v+f*o*o+m*o+n] += tmp2_[m*o*v*v+n*v*v+f*v+e];
                }
            }
        }
    }

    // - 1.00000 <e,i||a,b> t1(a,i) t2(b,f,m,n) 
    // + 1.00000 <f,i||a,b> t1(a,i) t2(b,e,m,n) 

#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int b = 0; b < v; b++) {
            double dum = 0.0;
            for (long int a = 0; a < v; a++) {
                for (long int i = 0; i < o; i++) {
                    dum += eri_aibc_[e*o*v*v+i*v*v+a*v+b] * t1_[a*o+i];
                }
            }
            tmp1_[e*v+b] = dum;
        }
    }
    F_DGEMM('n','n',o*o*v,v,v,1.0,t2_,o*o*v,tmp1_,v,0.0,tmp2_,o*o*v);
    C_DAXPY(o*o*v*v,-1.0,tmp2_,1,r2_,1);
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[e*o*o*v+f*o*o+m*o+n] += tmp2_[f*o*o*v+e*o*o+m*o+n];
                }
            }
        }
    }

    // - 1.00000 <e,i||a,b> t1(a,m) t2(b,f,n,i) 
    // + 1.00000 <f,i||a,b> t1(a,m) t2(b,e,n,i) 

    // v'(e,b,i,a) = <e,i||a,b>
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int i = 0; i < o; i++) {
            for (long int a = 0; a < v; a++) {
                for (long int b = 0; b < v; b++) {
                    tmp1_[e*o*v*v+b*o*v+i*v+a] = eri_aibc_[e*o*v*v+i*v*v+a*v+b];
                }
            }
        }
    }

    // I(m,e,b,i) = v'(e,b,i,a) t'(a,m)
    F_DGEMM('t','t',o*v*v,o,v,1.0,tmp1_,v,t1_,o,0.0,tmp2_,o*v*v);

    // t'(b,i,n,f) = t2(b,f,n,i)
#pragma omp parallel for schedule(static)
    for (long int b = 0; b < v; b++) {
        for (long int f = 0; f < v; f++) {
            for (long int n = 0; n < o; n++) {
                for (long int i = 0; i < o; i++) {
                    tmp1_[b*o*o*v+i*o*v+n*v+f] = t2_[b*o*o*v+f*o*o+n*o+i];
                }
            }
        }
    }

    // r'(m,e,n,f) = t'(b,i,n,f) I(m,e,b,i)
    F_DGEMM('n','n',o*v,o*v,o*v,1.0,tmp1_,o*v,tmp2_,o*v,0.0,tmp3_,o*v);

#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[e*o*o*v+f*o*o+m*o+n] -= tmp3_[m*o*v*v+e*o*v+n*v+f];
                    r2_[e*o*o*v+f*o*o+m*o+n] += tmp3_[m*o*v*v+f*o*v+n*v+e];
                }
            }
        }
    }

    // + 1.00000 <e,i||a,b> t1(a,n) ( t2(b,f,m,i) + t1(b,m) t1(f,i) )
    // - 1.00000 <f,i||a,b> t1(a,n) ( t2(b,e,m,i) + t1(b,m) t1(e,i) )

    // t'(b,i,m,f) = t2(b,f,m,i) + t1(b,m) t1(f,i)
#pragma omp parallel for schedule(static)
    for (long int b = 0; b < v; b++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int i = 0; i < o; i++) {
                    tmp1_[b*o*o*v+i*o*v+m*v+f] = t2_[b*o*o*v+f*o*o+m*o+i] + t1_[b*o+m] * t1_[f*o+i];
                }
            }
        }
    }

    // r(n,e,m,f) = t'(b,i,m,f) I(n,e,b,i) 
    F_DGEMM('n','n',o*v,o*v,o*v,1.0,tmp1_,o*v,tmp2_,o*v,0.0,tmp3_,o*v);

#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[e*o*o*v+f*o*o+m*o+n] += tmp3_[n*o*v*v+e*o*v+m*v+f];
                    r2_[e*o*o*v+f*o*o+m*o+n] -= tmp3_[n*o*v*v+f*o*v+m*v+e];
                }
            }
        }
    }

    //     + 0.25000 <i,j||a,b> t2(a,b,m,n) t2(e,f,i,j) 
    //     - 0.50000 <i,j||a,b> t1(a,n) t1(b,m) ( t2(e,f,i,j) + 2 t1(e,i) t1(f,j) )
    //     + 0.50000 <i,j||a,b> t1(e,i) t1(f,j) t2(a,b,m,n) 

    // t'(a,b,m,n) = t2(a,b,m,n) + 2 t1(a,m) t1(b,n)
    C_DCOPY(o*o*v*v,t2_,1,tmp1_,1);
#pragma omp parallel for schedule(static)
    for (long int a = 0; a < v; a++) {
        for (long int b = 0; b < v; b++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    tmp1_[a*o*o*v+b*o*o+m*o+n] += 2.0 * t1_[a*o+m] * t1_[b*o+n];
                }
            }
        }
    }

    // I(i,j,m,n) = t'(a,b,m,n) <i,j||a,b>
    F_DGEMM('n','n',o*o,o*o,v*v,1.0,tmp1_,o*o,eri_ijab_,v*v,0.0,tmp2_,o*o);

    // r(e,f,m,n) = I(i,j,m,n) t'(e,f,i,j)
    F_DGEMM('n','n',o*o,v*v,o*o,0.25,tmp2_,o*o,tmp1_,o*o,1.0,r2_,o*o);

    // - 0.50000 <i,j||a,b> t2(a,b,n,j) t2(e,f,m,i) 
    // + 0.50000 <i,j||a,b> t2(a,b,m,j) t2(e,f,n,i) 
    // - 1.00000 <i,j||a,b> t1(a,i) t1(b,n) t2(e,f,m,j) 
    // + 1.00000 <i,j||a,b> t1(a,i) t1(b,m) t2(e,f,n,j) 

    // t'(n,j,a,b) = 0.5 t2(a,b,n,j) + 2 t1(a,n) t1(b,j)
#pragma omp parallel for schedule(static)
    for (long int n = 0; n < o; n++) {
        for (long int j = 0; j < o; j++) {
            for (long int a = 0; a < v; a++) {
                for (long int b = 0; b < v; b++) {
                    tmp1_[n*o*v*v+j*v*v+a*v+b] = 0.5 * (t2_[a*o*o*v+b*o*o+n*o+j] + 2.0 * t1_[a*o+n] * t1_[b*o+j]);
                }
            }
        }
    }

    // I(i,n) = t'(n,j,a,b) <i,j||a,b>
    F_DGEMM('t','n',o,o,o*v*v,1.0,tmp1_,o*v*v,eri_ijab_,o*v*v,0.0,tmp2_,o);

    // r(e,f,m,n) = I(i,n) t2(e,f,m,i)
    F_DGEMM('n','n',o,o*v*v,o,1.0,tmp2_,o,t2_,o,0.0,tmp1_,o);

    C_DAXPY(o*o*v*v,-1.0,tmp1_,1,r2_,1);
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[e*o*o*v+f*o*o+m*o+n] += tmp1_[e*o*o*v+f*o*o+n*o+m];
                }
            }
        }
    }

    // - 0.50000 <i,j||a,b> t2(a,e,m,n) t2(b,f,i,j) 
    // + 0.50000 <i,j||a,b> t2(a,f,m,n) t2(b,e,i,j)
    // - 1.00000 <i,j||a,b> t1(a,i) t1(e,j) t2(b,f,m,n) 
    // + 1.00000 <i,j||a,b> t1(a,i) t1(f,j) t2(b,e,m,n) 
    
    // v'(a,b,i,j) = <i,j||a,b>
#pragma omp parallel for schedule(static)
    for (long int a = 0; a < v; a++) {
        for (long int b = 0; b < v; b++) {
            for (long int i = 0; i < o; i++) {
                for (long int j = 0; j < o; j++) {
                    tmp1_[a*o*o*v+b*o*o+i*o+j] = eri_ijab_[i*o*v*v+j*v*v+a*v+b];
                }
            }
        }
    }

    // t'(f,b,i,j) = 0.5 ( t2(b,f,i,j) + 2.0 t1(b,i) t1(f,j) )
#pragma omp parallel for schedule(static)
    for (long int f = 0; f < v; f++) {
        for (long int b = 0; b < v; b++) {
            for (long int i = 0; i < o; i++) {
                for (long int j = 0; j < o; j++) {
                    tmp2_[f*o*o*v+b*o*o+i*o+j] = 0.5 * ( t2_[b*o*o*v+f*o*o+i*o+j] + 2.0 * t1_[b*o+i] * t1_[f*o+j]);
                }
            }
        }
    }

    // I(f,a) = v'(a,b,i,j) t'(f,b,i,j)
    F_DGEMM('t','n',v,v,o*o*v,1.0,tmp1_,o*o*v,tmp2_,o*o*v,0.0,tmp3_,v);

    // r(f,e,m,n) = t2(a,e,m,n) I(f,a)
    F_DGEMM('n','n',o*o*v,v,v,1.0,t2_,o*o*v,tmp3_,v,0.0,tmp1_,o*o*v);

    C_DAXPY(o*o*v*v,1.0,tmp1_,1,r2_,1);
#pragma omp parallel for schedule(static)
    for (long int f = 0; f < v; f++) {
        for (long int e = 0; e < v; e++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[f*o*o*v+e*o*o+m*o+n] -= tmp1_[e*o*o*v+f*o*o+m*o+n];
                }   
            }   
        }   
    }

    // + 1.00000 <i,j||a,b> t2(a,e,n,j) t2(b,f,m,i)   (1st 2 can be unfolded long into 4 permutations)
    // - 1.00000 <i,j||a,b> t2(a,e,m,j) t2(b,f,n,i) 
    // + 1.00000 <i,j||a,b> t1(a,n) t1(e,j) t2(b,f,m,i) 
    // - 1.00000 <i,j||a,b> t1(a,n) t1(f,j) t2(b,e,m,i) 
    // - 1.00000 <i,j||a,b> t1(a,m) t1(e,j) t2(b,f,n,i) 
    // + 1.00000 <i,j||a,b> t1(a,m) t1(f,j) t2(b,e,n,i) 

    //I(i,b,j,a) = <i,j||a,b>
#pragma omp parallel for schedule(static)
    for (long int i = 0; i < o; i++) {
        for (long int j = 0; j < o; j++) {
            for (long int a = 0; a < v; a++) {
                for (long int b = 0; b < v; b++) {
                    tmp1_[i*o*v*v+b*o*v+j*v+a] = eri_ijab_[i*o*v*v+j*v*v+a*v+b];
                }
            }
        }
    }

    // t'(i,b,m,f) = t2(b,f,m,i)
#pragma omp parallel for schedule(static)
    for (long int b = 0; b < v; b++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int i = 0; i < o; i++) {
                    tmp2_[i*o*v*v+b*o*v+m*v+f] = t2_[b*o*o*v+f*o*o+m*o+i];
                }
            }
        }
    }

    //I''(m,f,j,a) = I(i,b,j,a) t'(i,b,m,f) 
    F_DGEMM('n','t',o*v,o*v,o*v,1.0,tmp1_,o*v,tmp2_,o*v,0.0,tmp3_,o*v);

    // t''(n,e,j,a) = t2(a,e,n,j) + 2 t1(a,n) t1(e,j)
#pragma omp parallel for schedule(static)
    for (long int n = 0; n < o; n++) {
        for (long int e = 0; e < v; e++) {
            for (long int j = 0; j < o; j++) {
                for (long int a = 0; a < v; a++) {
                    tmp1_[e*o*o*v+n*o*v+j*v+a] = t2_[a*o*o*v+e*o*o+n*o+j] + 2.0 * t1_[a*o+n] * t1_[e*o+j];
                }
            }
        }
    }

    //r'(e,n,m,f) = I''(m,f,j,a) t''(n,e,j,a)
    F_DGEMM('t','n',o*v,o*v,o*v,0.5,tmp3_,o*v,tmp1_,o*v,0.0,tmp2_,o*v);
#pragma omp parallel for schedule(static)
    for (long int e = 0; e < v; e++) {
        for (long int f = 0; f < v; f++) {
            for (long int m = 0; m < o; m++) {
                for (long int n = 0; n < o; n++) {
                    r2_[e*o*o*v+f*o*o+m*o+n] += tmp2_[e*o*o*v+n*o*v+m*v+f];
                    r2_[e*o*o*v+f*o*o+m*o+n] -= tmp2_[f*o*o*v+n*o*v+m*v+e];
                    r2_[e*o*o*v+f*o*o+m*o+n] -= tmp2_[e*o*o*v+m*o*v+n*v+f];
                    r2_[e*o*o*v+f*o*o+m*o+n] += tmp2_[f*o*o*v+m*o*v+n*v+e];
                }
            }
        }
    }

}

double PolaritonicUCCSD::update_amplitudes() {

    long int o = nalpha_ + nbeta_;
    long int v = 2 * nmo_ - o;

    // t2
    for (long int a = 0; a < v; a++) {
        double da = epsilon_[a+o];
        for (long int b = 0; b < v; b++) {
            double dab = da + epsilon_[b+o];
            for (long int i = 0; i < o; i++) {
                double dabi = dab - epsilon_[i];
                for (long int j = 0; j < o; j++) {
                    double dabij = dabi - epsilon_[j];
                    long int abij = a*o*o*v+b*o*o+i*o+j;
                    r2_[abij] /= -dabij;
                }
            }
        }
    }

    // t1
    for (long int a = 0; a < v; a++) {
        double da = epsilon_[a+o];
        for (long int i = 0; i < o; i++) {
            double dai = da - epsilon_[i];
            long int ai = a*o+i;
            r1_[ai] /= -dai;
        }
    }

    // diis
    diis->WriteVector(residual_);
    C_DAXPY(o*o*v*v+o*v,-1.0,tamps_,1,residual_,1);
    C_DAXPY(o*o*v*v+o*v,1.0,residual_,1,tamps_,1);
    diis->WriteErrorVector(residual_);
    diis->Extrapolate(tamps_);

    return C_DNRM2(o*o*v*v+o*v,residual_,1);

}

double PolaritonicUCCSD::correlation_energy() {

    long int o = nalpha_ + nbeta_;
    long int v = 2 * nmo_ - o;

    double ec = 0.0;

    // + 0.25000 <i,j||a,b> ( t2(a,b,i,j) + 2 * t1(a,i) t1(b,j) )
    for (long int a = 0; a < v; a++) {
        for (long int b = 0; b < v; b++) {
            for (long int i = 0; i < o; i++) {
                for (long int j = 0; j < o; j++) {
                    ec += 0.25 * eri_ijab_[i*o*v*v+j*v*v+a*v+b] * ( t2_[a*o*o*v+b*o*o+i*o+j] + 2.0 * t1_[a*o+i] * t1_[b*o+j] );
                }
            }
        }
    }

    return ec;

}

} // End namespaces
