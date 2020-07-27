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
#include <psi4/libpsi4util/process.h>
#include <psi4/libpsi4util/PsiOutStream.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/libqt/qt.h>

#include <math.h>

#include "jellium_integrals.h"
#include "jellium_scf_solver.h"

#include <misc/omp.h>
#include <misc/blas.h>
#include <misc/diis.h>

using namespace psi;

namespace hilbert{

Jellium_SCFSolver::Jellium_SCFSolver(Options & options)
    :options_(options){

    outfile->Printf("\n");
    outfile->Printf( "        ***************************************************\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    Jellium Hartree-Fock                         *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        ***************************************************\n");

    outfile->Printf("\n");

    // jellium integrals
    jelly_ = (std::shared_ptr<JelliumIntegrals>)(new JelliumIntegrals(options_));

    // paramters
    nso_       = options_.get_int("N_BASIS_FUNCTIONS"); 
    nelectron_ = options_.get_int("N_ELECTRONS");
    nsopi_     = jelly_->nsopi();
    doccpi_    = jelly_->doccpi();
    nirrep_    = jelly_->nirrep();

    // number of electrons
    if ( nelectron_ % 2 != 0 ) {
        throw PsiException("finite jellium scf currently requires an even number of electrons.",__FILE__,__LINE__);
    }

    // factor for box size ... hard coded for now to give <rho> = 1
    Lfac_      = pow(nelectron_,1.0/3.0) / M_PI;
    boxlength_ = M_PI / 2; 
    
    //since box is already pi a.u. long
    //double length_nm = options_.get_double("length");
    //Lfac_ = length_nm * 18.89725988 / M_PI;
    //boxlength_ = options_.get_double("LENGTH");

    // kinetic energy integrals
    T_ = jelly_->Ke;
    T_->scale(1.0/Lfac_/Lfac_);

    // potential energy integrals
    V_ = jelly_->NucAttrac;
    V_->scale(1.0/Lfac_);
    V_->scale(nelectron_); 

    // density matrix
    Da_ = (std::shared_ptr<Matrix>)(new Matrix(nirrep_,nsopi_,nsopi_));

    // so/mo transformation matrix
    Ca_ = (std::shared_ptr<Matrix>)(new Matrix(nirrep_,nsopi_,nsopi_));

    // fock matrix
    Fa_ = (std::shared_ptr<Matrix>)(new Matrix(nirrep_,nsopi_,nsopi_));

}

Jellium_SCFSolver::~Jellium_SCFSolver(){
}

double Jellium_SCFSolver::compute_energy(){

    // convergence parameters
    double e_convergence = options_.get_double("E_CONVERGENCE");
    double d_convergence = options_.get_double("D_CONVERGENCE");
    int    maxiter       = options_.get_int("MAXITER");
    bool   do_diis       = options_.get_bool("DIIS");

    // print some information about this computation
    outfile->Printf("\n");
    outfile->Printf("    ==> Jellium Hartree-Fock parameters <==\n");
    outfile->Printf("\n");
    outfile->Printf("    number of electrons:              %10i\n",nelectron_);
    outfile->Printf("    number of basis functions:        %10i\n",nso_);
    outfile->Printf("    maximum particle-in-a-box state:  %10i\n",jelly_->get_nmax());
    outfile->Printf("    e_convergence                     %10.2le\n",e_convergence);
    outfile->Printf("    d_convergence                     %10.2le\n",d_convergence);
    outfile->Printf("    maxiter                           %10i\n",maxiter);
    outfile->Printf("    diis?                             %10s\n",do_diis ? "yes" : "no");
    outfile->Printf("\n");

    // diis solver
    std::shared_ptr<DIIS> diis (new DIIS(nso_*nso_));
    
    // core hamiltonian
    std::shared_ptr<Matrix> h (new Matrix(T_));
    h->add(V_);

    // eigenvectors / eigenvalues of fock matrix
    std::shared_ptr<Vector> Feval (new Vector(nirrep_,nsopi_));

    // diagonalize core hamiltonian, get orbitals
    Fa_->diagonalize(Ca_,Feval);

    // build density matrix 
    for (int h = 0; h < nirrep_; h++) {
        C_DGEMM('n','t',nsopi_[h],nsopi_[h],doccpi_[h],1.0,
            &(Ca_->pointer(h)[0][0]),nsopi_[h],
            &(Ca_->pointer(h)[0][0]),nsopi_[h],0.0,
            &(Da_->pointer(h)[0][0]),nsopi_[h]);
    }

    double energy = Da_->vector_dot(h) + Da_->vector_dot(Fa_);
    outfile->Printf("    initial energy: %20.12lf\n",energy);
    outfile->Printf("\n");
    
    int iter = 0;
    
    outfile->Printf("    ");
    outfile->Printf("  iter");
    outfile->Printf("              energy");
    outfile->Printf("                |dE|");
    outfile->Printf("             |FD-DF|\n");
    
    double dele  = 0.0;
    double g_nrm = 0.0;
   
    // containers for J and K 
    std::shared_ptr<Matrix> Ja (new Matrix(nirrep_,nsopi_,nsopi_));
    std::shared_ptr<Matrix> Ka (new Matrix(nirrep_,nsopi_,nsopi_));

    do {

        build_J(Da_,Ja);
        build_K(Da_,Ka);

        // fock matrix
        Fa_->copy(Ja);
        Fa_->scale(2.0);
        Fa_->subtract(Ka);
        Fa_->scale(1.0/Lfac_);
        Fa_->add(h);

        // evaluate current energy, E = D(H+F)
        double new_energy = 0.0;
        new_energy += (nelectron_*nelectron_/2.0) * jelly_->selfval/Lfac_;
        new_energy += Da_->vector_dot(h);
        new_energy += Da_->vector_dot(Fa_);

        // orbital gradient
        std::shared_ptr<Matrix> orbital_gradient(new Matrix("[F,D]", nirrep_, nsopi_, nsopi_));
        orbital_gradient->gemm(false, false, 1.0, Fa_, Da_, 0.0);

        std::shared_ptr<Matrix> tmp(orbital_gradient->transpose());
        orbital_gradient->subtract(tmp);
        tmp.reset();

        g_nrm = 0.0;
        for (int h = 0; h < nirrep_; h++) {
            double dum = C_DNRM2(nsopi_[h]*nsopi_[h],orbital_gradient->pointer(h)[0],1);;
            g_nrm += dum * dum;
        }
        g_nrm = sqrt(g_nrm);

        // DIIS
        if (do_diis) {

            // diis wants data in vector
           
            std::shared_ptr<Vector> tmp_F(new Vector(nso_*nso_));
            std::shared_ptr<Vector> tmp_G(new Vector(nso_*nso_));

            int off = 0;
            for (int h = 0; h < nirrep_; h++) {
                C_DCOPY(nsopi_[h]*nsopi_[h],Fa_->pointer(h)[0],1,tmp_F->pointer() + off,1);
                C_DCOPY(nsopi_[h]*nsopi_[h],orbital_gradient->pointer(h)[0],1,tmp_G->pointer() + off,1);
                off += nsopi_[h]*nsopi_[h];
            }

            diis->WriteVector(tmp_F->pointer());
            diis->WriteErrorVector(tmp_G->pointer());
            diis->Extrapolate(tmp_F->pointer());

            off = 0;
            for (int h = 0; h < nirrep_; h++) {
                C_DCOPY(nsopi_[h]*nsopi_[h],tmp_F->pointer()+off,1,Fa_->pointer(h)[0],1);
                C_DCOPY(nsopi_[h]*nsopi_[h],tmp_G->pointer() + off,1,orbital_gradient->pointer(h)[0],1);
                off += nsopi_[h]*nsopi_[h];
            }

        }

        dele = fabs(new_energy - energy);

        // evaluate new orbitals
        std::shared_ptr<Vector> Feval (new Vector(nirrep_, nsopi_));
        Fa_->diagonalize(Ca_,Feval);

        // evaluate new density
        for (int h = 0; h < nirrep_; h++) {
            C_DGEMM('n','t',nsopi_[h],nsopi_[h],doccpi_[h],1.0,
                &(Ca_->pointer(h)[0][0]),nsopi_[h],
                &(Ca_->pointer(h)[0][0]),nsopi_[h],0.0,
                &(Da_->pointer(h)[0][0]),nsopi_[h]);
        }

        outfile->Printf("    %6i%20.12lf%20.12lf%20.12lf\n", iter, new_energy, dele, g_nrm);
        energy = new_energy;


        iter++;
        if (iter > maxiter) break;
       
    }while (dele > e_convergence || g_nrm > d_convergence);

    if (iter > maxiter) {
        throw PsiException("jellium scf did not converge.", __FILE__, __LINE__);
    }

    outfile->Printf("\n");
    outfile->Printf("      SCF iterations converged!\n");
    outfile->Printf("\n");

    //double fock_energy = D->vector_dot(K) / Lfac_;
    //outfile->Printf("      Fock energy:             %20.12lf\n", fock_energy);

    outfile->Printf("    * Jellium SCF total energy: %20.12lf\n", energy);
    outfile->Printf("\n");

    Process::environment.globals["CURRENT ENERGY"]    = energy;
    Process::environment.globals["JELLIUM SCF TOTAL ENERGY"] = energy;

    //CIS_slow();

    return energy;

}

void Jellium_SCFSolver::build_J(std::shared_ptr<Matrix> Da, std::shared_ptr<Matrix> Ja){

    for (int hp = 0; hp < nirrep_; hp++) {

        int offp = 0;
        double ** j_p = Ja->pointer(hp);
        for (int myh = 0; myh < hp; myh++) {
            offp += nsopi_[myh];
        }

        #pragma omp parallel for schedule(dynamic)
        for (int p = 0; p < nsopi_[hp]; p++) {
            int pp = p + offp;

            for (int q = p; q < nsopi_[hp]; q++) {
                int qq = q + offp;

                double myJ = 0.0;
                
                for (int hr = 0; hr < nirrep_; hr++) {
                    double ** d_p = Da->pointer(hr);
                    int offr = 0;
                    
                    for (int myh = 0; myh < hr; myh++) {
                        offr += nsopi_[myh];
                    }

                    for (int r = 0; r < nsopi_[hr]; r++) {
                        int rr = r + offr;
                        double dum = 0.0;
                        for (int s = r+1; s < nsopi_[hr]; s++) {
                            int ss = s + offr;
                            dum += d_p[r][s] * jelly_->ERI(pp,qq,rr,ss);
                        }
                        myJ += 2.0 * dum + d_p[r][r] * jelly_->ERI(pp,qq,rr,rr);
                    }
                }
                j_p[p][q] = myJ;
                j_p[q][p] = myJ;
            }
        }
    }

}

void Jellium_SCFSolver::build_K(std::shared_ptr<Matrix> Da, std::shared_ptr<Matrix> Ka){

    for (int hp = 0; hp < nirrep_; hp++) {

        int offp = 0;
        double ** k_p = Ka->pointer(hp);
        for (int myh = 0; myh < hp; myh++) {
            offp += nsopi_[myh];
        }

        #pragma omp parallel for schedule(dynamic)
        for (int p = 0; p < nsopi_[hp]; p++) {
            int pp = p + offp;

            for (int q = p; q < nsopi_[hp]; q++) {
                int qq = q + offp;

                double myK = 0.0;
                
                for (int hr = 0; hr < nirrep_; hr++) {
                    double ** d_p = Da->pointer(hr);
                    int offr = 0;
                    
                    for (int myh = 0; myh < hr; myh++) {
                        offr += nsopi_[myh];
                    }

                    // exchange
                    for (int r = 0; r < nsopi_[hr]; r++) {
                        int rr = r + offr;
                        for (int s = 0; s < nsopi_[hr]; s++) {
                            int ss = s + offr;
                            myK += d_p[r][s] * jelly_->ERI(pp,ss,rr,qq);
                        }
                    }
                }
                k_p[p][q] = myK;
                k_p[q][p] = myK;
            }
        }
    }
}

void Jellium_SCFSolver::CIS_slow() {

    int o = nelectron_ / 2;
    int v = nso_ - o;

    int * virpi = (int*)malloc(nirrep_*sizeof(int));
    int * ovpi  = (int*)malloc(nirrep_*sizeof(int));
    for (int h = 0; h < nirrep_; h++) {
        virpi[h] = nsopi_[h] - doccpi_[h];
        ovpi[h]  = 0;
    }
    for (int ho = 0; ho < nirrep_; ho++) {
        for (int hv = 0; hv < nirrep_; hv++) {
            int hov = ho ^ hv;
            ovpi[hov] += doccpi_[ho] * virpi[hv];
        }
    }

    /// transform fock matrix to mo basis
    Fa_->transform(Ca_);
    Fa_->print();

    double * Fij = (double*)malloc(o*o*sizeof(double));
    double * Fab = (double*)malloc(v*v*sizeof(double));

    memset((void*)Fij,'\0',o*o*sizeof(double));
    memset((void*)Fab,'\0',v*v*sizeof(double));

    int * symmetry = (int*)malloc(nso_*sizeof(int));

    int off = 0;
    for (int h = 0; h < nirrep_; h++) {
        double ** fp = Fa_->pointer(h);
        for (int i = 0; i < doccpi_[h]; i++) {
            int ii = i + off;
            for (int j = 0; j < doccpi_[h]; j++) {
                int jj = j + off;
                Fij[ii*o+jj] = fp[i][j];
            }
            symmetry[ii] = h;
        }
        off += doccpi_[h];
    }
    off = 0;
    for (int h = 0; h < nirrep_; h++) {
        double ** fp = Fa_->pointer(h);
        for (int a = 0; a < virpi[h]; a++) {
            int aa = a + off;
            for (int b = 0; b < virpi[h]; b++) {
                int bb = b + off;
                Fab[aa*v+bb] = fp[a+doccpi_[h]][b+doccpi_[h]];
            }
            symmetry[aa+o] = h;
        }
        off += virpi[h];
    }

    /// transform ERIs to mo basis

    std::shared_ptr<Matrix> iajb (new Matrix(nirrep_,ovpi,ovpi));
    std::shared_ptr<Matrix> ijab (new Matrix(nirrep_,ovpi,ovpi));

    int * ia = (int*)malloc(nirrep_*sizeof(int));
    int * jb = (int*)malloc(nirrep_*sizeof(int));

    // (ia|jb)
    memset((void*)ia,'\0',nirrep_*sizeof(int));
    for (int hi = 0; hi < nirrep_; hi++) {
        for (int ha = 0; ha < nirrep_; ha++) {
            int hia = hi ^ ha;
            if ( ovpi[hia] == 0 ) continue;
            for (int i = 0; i < doccpi_[hi]; i++) {
                double ** ci = Ca_->pointer(hi);
                for (int a = 0; a < virpi[ha]; a++) {
                    double ** ca = Ca_->pointer(ha);

                    double ** iajb_p = iajb->pointer(hia);
                    double ** ijab_p = ijab->pointer(hia);

                    memset((void*)jb,'\0',nirrep_*sizeof(int));
                    for (int hj = 0; hj < nirrep_; hj++) {
                        for (int hb = 0; hb < nirrep_; hb++) {
                            int hjb = hj ^ hb;
                            if ( hia != hjb ) continue;

                            for (int j = 0; j < doccpi_[hj]; j++) {
                                double ** cj = Ca_->pointer(hj);
                                for (int b = 0; b < virpi[hb]; b++) {
                                    double ** cb = Ca_->pointer(hb);

                                    // (ia|jb) = (mu nu | lambda sigma) c(i,mu) c(a,nu) c(j,lambda) c(b,sigma)
                                    double dum_iajb = 0.0;
                                    int offmu = 0;
                                    for (int mu = 0; mu < nsopi_[hi]; mu++) {
                                        int offnu = 0;
                                        for (int nu = 0; nu < nsopi_[ha]; nu++) {
                                            int offlambda = 0;
                                            for (int lambda = 0; lambda < nsopi_[hj]; lambda++) {
                                                int offsigma = 0;
                                                for (int sigma = 0; sigma < nsopi_[hb]; sigma++) {
                                                    dum_iajb += jelly_->ERI(mu+offmu,nu+offnu,lambda+offlambda,sigma+offsigma)
                                                              * ci[i            ][mu]
                                                              * ca[a+doccpi_[ha]][nu]
                                                              * cj[j            ][lambda]
                                                              * cb[b+doccpi_[hb]][sigma];
                                                }
                                                offsigma += nsopi_[hb];
                                            }
                                            offlambda += nsopi_[hj];
                                        }
                                        offnu += nsopi_[ha];
                                    }
                                    offmu += nsopi_[hi];
                                    iajb_p[ia[hia]][jb[hjb]] = dum_iajb;

                                    // (ij|ab) = (mu nu | lambda sigma) c(i,mu) c(j,nu) c(a,lambda) c(b,sigma)
                                    double dum_ijab = 0.0;
                                    offmu = 0;
                                    for (int mu = 0; mu < nsopi_[hi]; mu++) {
                                        int offnu = 0;
                                        for (int nu = 0; nu < nsopi_[hj]; nu++) {
                                            int offlambda = 0;
                                            for (int lambda = 0; lambda < nsopi_[ha]; lambda++) {
                                                int offsigma = 0;
                                                for (int sigma = 0; sigma < nsopi_[hb]; sigma++) {
                                                    dum_ijab += jelly_->ERI(mu+offmu,nu+offnu,lambda+offlambda,sigma+offsigma)
                                                              * ci[i            ][mu]
                                                              * cj[j            ][nu]
                                                              * ca[a+doccpi_[ha]][lambda]
                                                              * cb[b+doccpi_[hb]][sigma];
                                                    offsigma += nsopi_[hb];
                                                }
                                                offlambda += nsopi_[ha];
                                            }
                                            offnu += nsopi_[hj];
                                        }
                                        offmu += nsopi_[hi];
                                    }
                                    ijab_p[ia[hia]][jb[hjb]] = dum_ijab;

                                    jb[hjb]++;
                                }
                            }
                        }
                    }
                    ia[hia]++;
                }
            }
        }
    }

    // cis hamiltonian
    std::shared_ptr<Matrix> cis_ham (new Matrix(nirrep_,ovpi,ovpi));

    memset((void*)ia,'\0',nirrep_*sizeof(int));
    for (int hi = 0; hi < nirrep_; hi++) {
        for (int ha = 0; ha < nirrep_; ha++) {
            int hia = hi ^ ha;
            if ( ovpi[hia] == 0 ) continue;
            for (int i = 0; i < doccpi_[hi]; i++) {
                for (int a = 0; a < virpi[ha]; a++) {

                    double ** ham_p = cis_ham->pointer(hia);
                    double ** iajb_p = iajb->pointer(hia);
                    double ** ijab_p = ijab->pointer(hia);

                    memset((void*)jb,'\0',nirrep_*sizeof(int));
                    for (int hj = 0; hj < nirrep_; hj++) {
                        for (int hb = 0; hb < nirrep_; hb++) {
                            int hjb = hj ^ hb;
                            if ( hia != hjb ) continue;

                            for (int j = 0; j < doccpi_[hj]; j++) {
                                for (int b = 0; b < virpi[hb]; b++) {

                                    double dum = 2.0 * iajb_p[ia[hia]][jb[hjb]] - ijab_p[ia[hia]][jb[hjb]];
                                    if ( hi == hj && i == j ) {
                                        double ** fp = Fa_->pointer(ha);
                                        dum += fp[a+doccpi_[ha]][b+doccpi_[ha]];
                                    }
                                    if ( ha == hb && a == b ) {
                                        double ** fp = Fa_->pointer(hi);
                                        dum -= fp[i][j];
                                    }
                                    ham_p[ia[hia]][jb[hjb]] = dum;
                                    
                                    jb[hjb]++;
                                }
                            }
                        }
                    }
                    ia[hia]++;
                }
            }
        }
    }

    std::shared_ptr<Matrix> cis_eigvec (new Matrix(cis_ham));
    std::shared_ptr<Vector> cis_eigval (new Vector(nirrep_,ovpi));
    cis_ham->diagonalize(cis_eigvec,cis_eigval,ascending);
    cis_eigval->print();

    free(virpi);
    free(ovpi);
}



} // End namespaces



