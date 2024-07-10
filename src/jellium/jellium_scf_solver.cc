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

//#include <algorithm>

#include <psi4/psi4-dec.h>
#include <psi4/physconst.h>
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
#include <misc/davidson_solver.h>

using namespace psi;

namespace hilbert{

static void evaluate_sigma(size_t N, size_t maxdim,double **sigma, double **b, void * data) {
    Jellium_SCFSolver* jellium = reinterpret_cast<Jellium_SCFSolver*>(data);
    jellium->CIS_evaluate_sigma(N,maxdim,b,sigma);
}
static double hamiltonian_element(size_t i, size_t j, void * data) {
    Jellium_SCFSolver* jellium = reinterpret_cast<Jellium_SCFSolver*>(data);
    return jellium->evaluate_hamiltonian_element(i,j);
}

Jellium_SCFSolver::Jellium_SCFSolver(Options & options)
    :options_(options){

    outfile->Printf("\n");
    outfile->Printf( "        ***************************************************\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    Finite Jellium Hartree-Fock                  *\n");
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
    if ( !options_["JELLIUM_BOX_LENGTH"].has_changed() ) {
        boxlength_ = pow(nelectron_,1.0/3.0);
        Lfac_ = boxlength_ / M_PI;
    }else {
        //since box is already pi a.u. long
        double length_nm = options_.get_double("JELLIUM_BOX_LENGTH");
        boxlength_ = length_nm * 10.0 / pc_bohr2angstroms;
        Lfac_ = boxlength_ / M_PI;
    }

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

    if ( options_["DOCC"].has_changed() ) {
        if ( options_["DOCC"].size() != nirrep_ ) {
            throw PSIEXCEPTION("Input DOCC array has the wrong dimensions");
        }
        for (int h = 0; h < nirrep_; h++) {
            doccpi_[h] = options_["DOCC"][h].to_integer();
        }
    }
    print_occupations();

    // diis solver
    std::shared_ptr<DIIS> diis (new DIIS(nso_*nso_));
    
    // core hamiltonian
    std::shared_ptr<Matrix> h (new Matrix(T_));
    h->add(V_);

    // eigenvectors / eigenvalues of fock matrix
    Dimension nsopi(nirrep_);
    for (int h = 0; h < nirrep_; h++){
        nsopi[h] = nsopi_[h];
    }
    std::shared_ptr<Vector> epsilon_a = std::make_shared<Vector>(nsopi);

    // diagonalize core hamiltonian, get orbitals
    Fa_->diagonalize(Ca_,epsilon_a);

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
        Fa_->diagonalize(Ca_,epsilon_a);

        // update occupations
        update_occupations(epsilon_a);

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

    outfile->Printf("    * Jellium SCF total energy: %20.12lf\n", energy);
    outfile->Printf("\n");

    print_occupations();

    Process::environment.globals["CURRENT ENERGY"]    = energy;
    Process::environment.globals["JELLIUM SCF TOTAL ENERGY"] = energy;

    //  O(n^5) scaling integral transformation, 2 N^4 + o^2 v^2 storage
    CIS_in_core();

    //  O(n^6) scaling integral transformation, o v N^2 + o^2 v^2 storage
    //CIS_in_core_n6();

    //  O(n^8) scaling integral transformation, o^2 v^2 storage
    //CIS_slow();

    // O(N^4 * m) scaling, m = number of desired roots
    //CIS_direct();


    return energy;

}

void Jellium_SCFSolver::update_occupations(std::shared_ptr<Vector> epsilon_a) {

    if ( options_["DOCC"].has_changed() ) {
        return;
    }

    int * old_docc = (int*)malloc(nirrep_*sizeof(int));
    for (int h = 0; h < nirrep_; ++h) {
        old_docc[h] = doccpi_[h];
    }

    std::vector<std::pair<double, int> > pairs_a;
    for (int h = 0; h < epsilon_a->nirrep(); ++h) {
        for (int i = 0; i < epsilon_a->dimpi()[h]; ++i) {
            pairs_a.push_back(std::make_pair(epsilon_a->get(h, i), h));
        }
    }
    sort(pairs_a.begin(), pairs_a.end());

    memset(doccpi_, 0, sizeof(int) * epsilon_a->nirrep());
    for (int i = 0; i < nelectron_ / 2; ++i) doccpi_[pairs_a[i].second]++;

    bool occ_changed = false;
    for (int h = 0; h < nirrep_; ++h) {
        if (old_docc[h] != doccpi_[h]) {
            occ_changed = true;
            break;
        }
    }

    if (occ_changed) {
        print_occupations();
    }

    free(old_docc);
}

void Jellium_SCFSolver::print_occupations() {
    outfile->Printf("\n");
    outfile->Printf("    Occupation by irrep:\n");
    outfile->Printf("\n");
    outfile->Printf("           ");
    for (int h = 0; h < nirrep_; ++h) outfile->Printf(" %4s ", jelly_->labels[h].c_str());
    outfile->Printf("\n");
    outfile->Printf("    DOCC [ ");
    for (int h = 0; h < nirrep_ - 1; ++h) outfile->Printf(" %4d,", doccpi_[h]);
    outfile->Printf(" %4d ]\n", doccpi_[nirrep_ - 1]);
    outfile->Printf("\n");
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


void Jellium_SCFSolver::CIS_in_core() {

    outfile->Printf("\n");
    outfile->Printf( "        ***************************************************\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    Finite Jellium CIS                           *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        ***************************************************\n");

    outfile->Printf("\n");
    int o = nelectron_ / 2;
    int v = nso_ - o;

    Dimension virpi(nirrep_);
    Dimension ovpi(nirrep_);
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
    std::shared_ptr<Matrix> F (new Matrix(Fa_));
    F->transform(Ca_);

    // dipole integrals
    std::shared_ptr<Vector> dipole_x = std::make_shared<Vector>(ovpi);
    std::shared_ptr<Vector> dipole_y = std::make_shared<Vector>(ovpi);
    std::shared_ptr<Vector> dipole_z = std::make_shared<Vector>(ovpi);

    int * symmetry = (int*)malloc(nso_*sizeof(int));
    for (int h = 0; h < nirrep_; h++) {
        int off = 0; 
        for (int myh = 0; myh < h; myh++) {
            off += nsopi_[myh];
        }
        for (int i = 0; i < nsopi_[h]; i++) {
            symmetry[i+off] = h;
        }
    }

    /// transform ERIs to mo basis

    outfile->Printf("    transform (ia|jb), (ij|ab)...."); fflush(stdout);
    std::shared_ptr<Matrix> cis_ham (new Matrix(nirrep_,ovpi,ovpi));
    for (int h = 0; h < nirrep_; h++) {

        int num_roots = 10;
        if ( options_["ROOTS_PER_IRREP"].has_changed() ) {
            num_roots = options_["ROOTS_PER_IRREP"][h].to_integer();
        }
        if ( num_roots == 0 ) continue;

        if ( ovpi[h] == 0 ) continue;

        // list of CIS transitions in this irrep
        cis_transition_list_.clear();
        for (int hi = 0; hi < nirrep_; hi++) {
            for (int ha = 0; ha < nirrep_; ha++) {
                int hia = hi ^ ha;
                if ( hia != h ) continue;
                for (int i = 0; i < doccpi_[hi]; i++) {
                    for (int a = 0; a < virpi[ha]; a++) {

                        cis_transition me;
                        me.i   = i;
                        me.a   = a;
                        me.hi  = hi;
                        me.ha  = ha;

                        cis_transition_list_.push_back(me);
                    }
                }
            }
        }

        // dipole integrals
        double * dipole_x_p = dipole_x->pointer(h);
        double * dipole_y_p = dipole_y->pointer(h);
        double * dipole_z_p = dipole_z->pointer(h);
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** ci = Ca_->pointer(hi);
            double ** ca = Ca_->pointer(ha);

            double dum_x = 0.0;
            double dum_y = 0.0;
            double dum_z = 0.0;
            for (int mu = 0; mu < nsopi_[hi]; mu++) {
                for (int nu = 0; nu < nsopi_[ha]; nu++) {

                    double dip_x = jelly_->dipole_x(mu+i_off,nu+a_off,boxlength_);
                    double dip_y = jelly_->dipole_y(mu+i_off,nu+a_off,boxlength_);
                    double dip_z = jelly_->dipole_z(mu+i_off,nu+a_off,boxlength_);
                    dum_x += dip_x
                           * ci[mu][i            ]
                           * ca[nu][a+doccpi_[ha]];
                    dum_y += dip_y
                           * ci[mu][i            ]
                           * ca[nu][a+doccpi_[ha]];
                    dum_z += dip_z
                           * ci[mu][i            ]
                           * ca[nu][a+doccpi_[ha]];
                }
                dipole_x_p[ia] = dum_x;
                dipole_y_p[ia] = dum_y;
                dipole_z_p[ia] = dum_z;
            }
        }

        double ** ham_p = cis_ham->pointer(h);

        // todo: make these smaller
        double * tmp_1 = (double*)malloc(nso_*nso_*nso_*nso_*sizeof(double));
        double * tmp_2 = (double*)malloc(nso_*nso_*nso_*nso_*sizeof(double));
        memset((void*)tmp_1,'\0',nso_*nso_*nso_*nso_*sizeof(double));
        memset((void*)tmp_2,'\0',nso_*nso_*nso_*nso_*sizeof(double));
    
        // (ia|jb)

        // (i nu| lambda sigma) = c(mu,i) (mu nu | lambda sigma) first partial transformation
        for (int hi = 0; hi < nirrep_; hi++) {
            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }
            double ** ci = Ca_->pointer(hi);
            for (int i = 0; i < doccpi_[hi]; i++) {
                for (int h_nu = 0; h_nu < nirrep_; h_nu++) {
                    int hi_nu = hi ^ h_nu;
                    if ( hi_nu != h ) continue;
                    int nu_off = 0;
                    for (int myh = 0; myh < h_nu; myh++) {
                        nu_off += nsopi_[myh];
                    }
                    #pragma omp parallel for schedule(dynamic)
                    for (int nu = 0; nu < nsopi_[h_nu]; nu++) {

                        for (int h_lambda = 0; h_lambda < nirrep_; h_lambda++) {
                            int lambda_off = 0;
                            for (int myh = 0; myh < h_lambda; myh++) {
                                lambda_off += nsopi_[myh];
                            }
                            for (int h_sigma = 0; h_sigma < nirrep_; h_sigma++) {
                                int h_lambda_sigma = h_lambda ^ h_sigma;
                                if ( h_lambda_sigma != h ) continue;
                                int sigma_off = 0;
                                for (int myh = 0; myh < h_sigma; myh++) {
                                    sigma_off += nsopi_[myh];
                                }
                                for (int lambda = 0; lambda < nsopi_[h_lambda]; lambda++) {
                                    for (int sigma = 0; sigma < nsopi_[h_sigma]; sigma++) {

                                        // sum over mu
                                        double dum = 0.0;
                                        for (int mu = 0; mu < nsopi_[hi]; mu++) {
                                            dum += jelly_->ERI(mu+i_off,nu+nu_off,lambda+lambda_off,sigma+sigma_off) * ci[mu][i];
                                        }
                                        tmp_1[(i+i_off)*nso_*nso_*nso_ + (nu+nu_off)*nso_*nso_ + (lambda+lambda_off)*nso_ + (sigma+sigma_off)] = dum;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // (i nu| j sigma) = c(lambda,j) (i nu | lambda sigma) second partial transformation
        for (int hi = 0; hi < nirrep_; hi++) {
            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }
            double ** ci = Ca_->pointer(hi);
            for (int i = 0; i < doccpi_[hi]; i++) {
                for (int h_nu = 0; h_nu < nirrep_; h_nu++) {
                    int hi_nu = hi ^ h_nu;
                    if ( hi_nu != h ) continue;
                    int nu_off = 0;
                    for (int myh = 0; myh < h_nu; myh++) {
                        nu_off += nsopi_[myh];
                    }
                    #pragma omp parallel for schedule(dynamic)
                    for (int nu = 0; nu < nsopi_[h_nu]; nu++) {

                        for (int hj = 0; hj < nirrep_; hj++) {
                            int j_off = 0;
                            for (int myh = 0; myh < hj; myh++) {
                                j_off += nsopi_[myh];
                            }
                            double ** cj = Ca_->pointer(hj);
                            for (int j = 0; j < doccpi_[hj]; j++) {

                                for (int h_sigma = 0; h_sigma < nirrep_; h_sigma++) {
                                    int h_j_sigma = hj ^ h_sigma;
                                    if ( h_j_sigma != h ) continue;
                                    int sigma_off = 0;
                                    for (int myh = 0; myh < h_sigma; myh++) {
                                        sigma_off += nsopi_[myh];
                                    }
                                    for (int sigma = 0; sigma < nsopi_[h_sigma]; sigma++) {
                                        // sum over lambda
                                        double dum = 0.0;
                                        for (int lambda = 0; lambda < nsopi_[hj]; lambda++) {
                                            dum += tmp_1[(i+i_off)*nso_*nso_*nso_ + (nu+nu_off)*nso_*nso_ + (lambda+j_off)*nso_ + (sigma+sigma_off)] * cj[lambda][j];
                                        }
                                        tmp_2[(i+i_off)*nso_*nso_*nso_ + (nu+nu_off)*nso_*nso_ + (j+j_off)*nso_ + (sigma+sigma_off)] = dum;
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // (i a| j sigma) = c(nu,a) (i nu | j sigma) third partial transformation
        #pragma omp parallel for schedule(dynamic)
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;
            int hia = hi ^ ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** ca = Ca_->pointer(ha);

            for (int hj = 0; hj < nirrep_; hj++) {
                int j_off = 0;
                for (int myh = 0; myh < hj; myh++) {
                    j_off += nsopi_[myh];
                }
                double ** cj = Ca_->pointer(hj);
                for (int j = 0; j < doccpi_[hj]; j++) {

                    for (int h_sigma = 0; h_sigma < nirrep_; h_sigma++) {
                        int h_j_sigma = hj ^ h_sigma;
                        if ( h_j_sigma != h ) continue;
                        int sigma_off = 0;
                        for (int myh = 0; myh < h_sigma; myh++) {
                            sigma_off += nsopi_[myh];
                        }
                        for (int sigma = 0; sigma < nsopi_[h_sigma]; sigma++) {

                            // sum over nu
                            double dum = 0.0;
                            for (int nu = 0; nu < nsopi_[ha]; nu++) {
                                dum += tmp_2[(i+i_off)*nso_*nso_*nso_ + (nu+a_off)*nso_*nso_ + (j+j_off)*nso_ + (sigma+sigma_off)] * ca[nu][a+doccpi_[ha]];
                            }
                            tmp_1[(i+i_off)*nso_*nso_*nso_ + (a+a_off)*nso_*nso_ + (j+j_off)*nso_ + (sigma+sigma_off)] = dum;

                        }
                    }
                }
            }
        }
        // (i a| j b) = c(sigma,b) (i a | j sigma) fourth partial transformation
        #pragma omp parallel for schedule(dynamic)
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            for (int jb = 0; jb < cis_transition_list_.size(); jb++) {

                int j   = cis_transition_list_[jb].i;
                int b   = cis_transition_list_[jb].a;
                int hj  = cis_transition_list_[jb].hi;
                int hb  = cis_transition_list_[jb].ha;

                int j_off = 0;
                for (int myh = 0; myh < hj; myh++) {
                    j_off += nsopi_[myh];
                }

                int b_off = 0;
                for (int myh = 0; myh < hb; myh++) {
                    b_off += nsopi_[myh];
                }

                double ** cb = Ca_->pointer(hb);
                // sum over sigma
                double dum = 0.0;
                for (int sigma = 0; sigma < nsopi_[hb]; sigma++) {
                    dum += tmp_1[(i+i_off)*nso_*nso_*nso_ + (a+a_off)*nso_*nso_ + (j+j_off)*nso_ + (sigma+b_off)] * cb[sigma][b+doccpi_[hb]];
                }
                ham_p[ia][jb] = 2.0 * dum / Lfac_;
            }
        }

        // (ij|ab)

        // (i nu| lambda sigma) = c(mu,i) (mu nu | lambda sigma) first partial transformation
        for (int hi = 0; hi < nirrep_; hi++) {
            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }
            double ** ci = Ca_->pointer(hi);
            for (int i = 0; i < doccpi_[hi]; i++) {

                for (int h_lambda = 0; h_lambda < nirrep_; h_lambda++) {
                    int hi_lambda = hi ^ h_lambda;
                    if ( hi_lambda != h ) continue;
                    int lambda_off = 0;
                    for (int myh = 0; myh < h_lambda; myh++) {
                        lambda_off += nsopi_[myh];
                    }
                    #pragma omp parallel for schedule(dynamic)
                    for (int lambda = 0; lambda < nsopi_[h_lambda]; lambda++) {

                        for (int h_nu = 0; h_nu < nirrep_; h_nu++) {
                            int nu_off = 0;
                            for (int myh = 0; myh < h_nu; myh++) {
                                nu_off += nsopi_[myh];
                            }
                            for (int h_sigma = 0; h_sigma < nirrep_; h_sigma++) {
                                int h_nu_sigma = h_nu ^ h_sigma;
                                if ( h_nu_sigma != h ) continue;
                                int sigma_off = 0;
                                for (int myh = 0; myh < h_sigma; myh++) {
                                    sigma_off += nsopi_[myh];
                                }
                                for (int nu = 0; nu < nsopi_[h_nu]; nu++) {
                                    for (int sigma = 0; sigma < nsopi_[h_sigma]; sigma++) {

                                        // sum over mu
                                        double dum = 0.0;
                                        for (int mu = 0; mu < nsopi_[hi]; mu++) {
                                            dum += jelly_->ERI(mu+i_off,nu+nu_off,lambda+lambda_off,sigma+sigma_off) * ci[mu][i];
                                        }
                                        tmp_1[(i+i_off)*nso_*nso_*nso_ + (nu+nu_off)*nso_*nso_ + (lambda+lambda_off)*nso_ + (sigma+sigma_off)] = dum;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // (i j| lambda sigma) = c(nu,j) (i nu | lambda sigma) second partial transformation
        for (int hi = 0; hi < nirrep_; hi++) {
            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }
            for (int i = 0; i < doccpi_[hi]; i++) {

                for (int h_lambda = 0; h_lambda < nirrep_; h_lambda++) {
                    int hi_lambda = hi ^ h_lambda;
                    if ( hi_lambda != h ) continue;
                    int lambda_off = 0;
                    for (int myh = 0; myh < h_lambda; myh++) {
                        lambda_off += nsopi_[myh];
                    }
                    #pragma omp parallel for schedule(dynamic)
                    for (int lambda = 0; lambda < nsopi_[h_lambda]; lambda++) {

                        for (int hj = 0; hj < nirrep_; hj++) {
                            int j_off = 0;
                            for (int myh = 0; myh < hj; myh++) {
                                j_off += nsopi_[myh];
                            }
                            double ** cj = Ca_->pointer(hj);
                            for (int j = 0; j < doccpi_[hj]; j++) {

                                for (int h_sigma = 0; h_sigma < nirrep_; h_sigma++) {
                                    int h_j_sigma = hj ^ h_sigma;
                                    if ( h_j_sigma != h ) continue;
                                    int sigma_off = 0;
                                    for (int myh = 0; myh < h_sigma; myh++) {
                                        sigma_off += nsopi_[myh];
                                    }
                                    for (int sigma = 0; sigma < nsopi_[h_sigma]; sigma++) {

                                        // sum over nu
                                        double dum = 0.0;
                                        for (int nu = 0; nu < nsopi_[hj]; nu++) {
                                            dum += tmp_1[(i+i_off)*nso_*nso_*nso_ + (nu+j_off)*nso_*nso_ + (lambda+lambda_off)*nso_ + (sigma+sigma_off)] * cj[nu][j];
                                        }
                                        tmp_2[(i+i_off)*nso_*nso_*nso_ + (j+j_off)*nso_*nso_ + (lambda+lambda_off)*nso_ + (sigma+sigma_off)] = dum;
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // (i j| a sigma) = c(lambda,a) (i j | lambda sigma) third partial transformation
        #pragma omp parallel for schedule(dynamic)
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;
            int hia = hi ^ ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** ca = Ca_->pointer(ha);

            for (int hj = 0; hj < nirrep_; hj++) {
                int j_off = 0;
                for (int myh = 0; myh < hj; myh++) {
                    j_off += nsopi_[myh];
                }
                double ** cj = Ca_->pointer(hj);
                for (int j = 0; j < doccpi_[hj]; j++) {

                    for (int h_sigma = 0; h_sigma < nirrep_; h_sigma++) {
                        int h_j_sigma = hj ^ h_sigma;
                        if ( h_j_sigma != h ) continue;
                        int sigma_off = 0;
                        for (int myh = 0; myh < h_sigma; myh++) {
                            sigma_off += nsopi_[myh];
                        }
                        for (int sigma = 0; sigma < nsopi_[h_sigma]; sigma++) {

                            // sum over lambda
                            double dum = 0.0;
                            for (int lambda = 0; lambda < nsopi_[ha]; lambda++) {
                                dum += tmp_2[(i+i_off)*nso_*nso_*nso_ + (j+j_off)*nso_*nso_ + (lambda+a_off)*nso_ + (sigma+sigma_off)] * ca[lambda][a+doccpi_[ha]];
                            }
                            tmp_1[(i+i_off)*nso_*nso_*nso_ + (j+j_off)*nso_*nso_ + (a+a_off)*nso_ + (sigma+sigma_off)] = dum;

                        }
                    }
                }
            }
        }
        // (i j| a b) = c(sigma,b) (i j | a b) fourth partial transformation
        #pragma omp parallel for schedule(dynamic)
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            for (int jb = 0; jb < cis_transition_list_.size(); jb++) {

                int j   = cis_transition_list_[jb].i;
                int b   = cis_transition_list_[jb].a;
                int hj  = cis_transition_list_[jb].hi;
                int hb  = cis_transition_list_[jb].ha;

                int j_off = 0;
                for (int myh = 0; myh < hj; myh++) {
                    j_off += nsopi_[myh];
                }

                int b_off = 0;
                for (int myh = 0; myh < hb; myh++) {
                    b_off += nsopi_[myh];
                }

                double ** cb = Ca_->pointer(hb);
                // sum over sigma
                double dum = 0.0;
                for (int sigma = 0; sigma < nsopi_[hb]; sigma++) {
                    dum += tmp_1[(i+i_off)*nso_*nso_*nso_ + (j+j_off)*nso_*nso_ + (a+a_off)*nso_ + (sigma+b_off)] * cb[sigma][b+doccpi_[hb]];
                }
                ham_p[ia][jb] -= dum / Lfac_;
            }
        }
        free(tmp_1);
        free(tmp_2);

        // remaining diagonal parts of Hamiltonian
        #pragma omp parallel for schedule(dynamic)
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** fa = F->pointer(ha);
            double ** fi = F->pointer(hi);
            ham_p[ia][ia] += fa[a+doccpi_[ha]][a+doccpi_[ha]] - fi[i][i];
        }
    }
    outfile->Printf("done.\n");

    free(symmetry);

    outfile->Printf("    diagonalize CIS Hamiltonian..."); fflush(stdout);
    std::shared_ptr<Matrix> cis_eigvec (new Matrix(cis_ham));
    std::shared_ptr<Vector> cis_eigval = std::make_shared<Vector>(ovpi);
    cis_ham->diagonalize(cis_eigvec,cis_eigval,ascending);
    outfile->Printf("done.\n");

    outfile->Printf("\n");
    outfile->Printf("    ==> CIS excitaiton energies <=="); fflush(stdout);
    outfile->Printf("\n");

    for (int h = 0; h < nirrep_; h++) {

        int num_roots = 10;
        if ( options_["ROOTS_PER_IRREP"].has_changed() ) {
            num_roots = options_["ROOTS_PER_IRREP"][h].to_integer();
        }
        if ( num_roots == 0 ) continue;

        outfile->Printf("\n");
        outfile->Printf("    Irrep: %10s\n",jelly_->labels[h].c_str());
        outfile->Printf("\n");
        outfile->Printf("    state");
        outfile->Printf("                w(Eh)");
        outfile->Printf("                w(eV)");
        outfile->Printf("                    f");
        outfile->Printf("\n");
        double * eigval_p = cis_eigval->pointer(h);
        for (int I = 0; I < ovpi[h]; I++) {

            // evaluate oscillator strengths
            double tdp_x = sqrt(2.0) * C_DDOT(cis_transition_list_.size(),&(cis_eigvec->pointer(h)[0][I]),cis_transition_list_.size(),dipole_x->pointer(h),1);
            double tdp_y = sqrt(2.0) * C_DDOT(cis_transition_list_.size(),&(cis_eigvec->pointer(h)[0][I]),cis_transition_list_.size(),dipole_y->pointer(h),1);
            double tdp_z = sqrt(2.0) * C_DDOT(cis_transition_list_.size(),&(cis_eigvec->pointer(h)[0][I]),cis_transition_list_.size(),dipole_z->pointer(h),1);
            double f = 2.0 / 3.0 * eigval_p[I] * ( tdp_x * tdp_x + tdp_y * tdp_y + tdp_z * tdp_z );

            outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",I,eigval_p[I],eigval_p[I] * pc_hartree2ev, f);
        }
        outfile->Printf("\n");

    }

    // check sum
    double check = 0.0;
    for (int h = 0; h < nirrep_; h++) {
        check += C_DNRM2(ovpi[h],cis_eigval->pointer(h),1);
    }
    outfile->Printf("\n");
    outfile->Printf("    ||eps|| = %20.12lf\n",check);
    outfile->Printf("\n");

    //free(virpi);
    //free(ovpi);
}
void Jellium_SCFSolver::CIS_in_core_n6() {

    outfile->Printf("\n");
    outfile->Printf( "        ***************************************************\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    Finite Jellium CIS                           *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        ***************************************************\n");

    outfile->Printf("\n");
    int o = nelectron_ / 2;
    int v = nso_ - o;

    Dimension virpi(nirrep_);
    Dimension ovpi(nirrep_);
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
    std::shared_ptr<Matrix> F (new Matrix(Fa_));
    F->transform(Ca_);

/*
    bool * skip = (bool*)malloc(nso_*sizeof(bool));
    memset((void*)skip,'\0',nso_*sizeof(bool));
    for (int i = 0; i < nso_; i++) {
        
        double min = 9e9;
        int minj   = 0;
        int minh   = 0;
        int minoff = 0;

        for (int hj = 0; hj < nirrep_; hj++) {
            double ** fp = F->pointer(hj);
            for (int j = 0; j < nsopi_[hj]; j++) {

                int j_off = 0;
                for (int myh = 0; myh < hj; myh++) {
                    j_off += nsopi_[myh];
                }

                if ( skip[j+j_off] ) continue;

                if ( fp[j][j] < min ) {
                    min  = fp[j][j];
                    minh = hj;
                    minj = j;
                    minoff = j_off;
                }
            }
        }
        printf("%5i %20.12lf\n",i,min);//F->pointer(minh)[minj][minj]);
        skip[minj+minoff] = true;
    }
    exit(0);
*/

    // dipole integrals
    std::shared_ptr<Vector> dipole_x = std::make_shared<Vector>(ovpi);
    std::shared_ptr<Vector> dipole_y = std::make_shared<Vector>(ovpi);
    std::shared_ptr<Vector> dipole_z = std::make_shared<Vector>(ovpi);

    /// transform ERIs to mo basis

    outfile->Printf("    transform (ia|jb), (ij|ab)...."); fflush(stdout);
    std::shared_ptr<Matrix> cis_ham (new Matrix(nirrep_,ovpi,ovpi));
    for (int h = 0; h < nirrep_; h++) {

        int num_roots = 10;
        if ( options_["ROOTS_PER_IRREP"].has_changed() ) {
            num_roots = options_["ROOTS_PER_IRREP"][h].to_integer();
        }
        if ( num_roots == 0 ) continue;

        if ( ovpi[h] == 0 ) continue;

        // list of CIS transitions in this irrep
        cis_transition_list_.clear();
        for (int hi = 0; hi < nirrep_; hi++) {
            for (int ha = 0; ha < nirrep_; ha++) {
                int hia = hi ^ ha;
                if ( hia != h ) continue;
                for (int i = 0; i < doccpi_[hi]; i++) {
                    for (int a = 0; a < virpi[ha]; a++) {

                        cis_transition me;
                        me.i   = i;
                        me.a   = a;
                        me.hi  = hi;
                        me.ha  = ha;

                        cis_transition_list_.push_back(me);
                    }
                }
            }
        }

        // dipole integrals
        double * dipole_x_p = dipole_x->pointer(h);
        double * dipole_y_p = dipole_y->pointer(h);
        double * dipole_z_p = dipole_z->pointer(h);
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** ci = Ca_->pointer(hi);
            double ** ca = Ca_->pointer(ha);

            double dum_x = 0.0;
            double dum_y = 0.0;
            double dum_z = 0.0;
            for (int mu = 0; mu < nsopi_[hi]; mu++) {
                for (int nu = 0; nu < nsopi_[ha]; nu++) {

                    double dip_x = jelly_->dipole_x(mu+i_off,nu+a_off,boxlength_);
                    double dip_y = jelly_->dipole_y(mu+i_off,nu+a_off,boxlength_);
                    double dip_z = jelly_->dipole_z(mu+i_off,nu+a_off,boxlength_);
                    dum_x += dip_x
                           * ci[mu][i            ]
                           * ca[nu][a+doccpi_[ha]];
                    dum_y += dip_y
                           * ci[mu][i            ]
                           * ca[nu][a+doccpi_[ha]];
                    dum_z += dip_z
                           * ci[mu][i            ]
                           * ca[nu][a+doccpi_[ha]];
                }
                dipole_x_p[ia] = dum_x;
                dipole_y_p[ia] = dum_y;
                dipole_z_p[ia] = dum_z;
            }
        }

        double ** ham_p = cis_ham->pointer(h);

        std::shared_ptr<Matrix> tmp (new Matrix(ovpi[h],nso_*nso_));
        double ** tmp_p = tmp->pointer();

        // (ia|jb) first half-transformation
        #pragma omp parallel for schedule(dynamic)
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;
            int hia = hi ^ ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** ci = Ca_->pointer(hi);
            double ** ca = Ca_->pointer(ha);

            for (int h_lambda = 0; h_lambda < nirrep_; h_lambda++) {
                int lambda_off = 0;
                for (int myh = 0; myh < h_lambda; myh++) {
                    lambda_off += nsopi_[myh];
                }
                for (int h_sigma = 0; h_sigma < nirrep_; h_sigma++) {
                    int sigma_off = 0;
                    for (int myh = 0; myh < h_sigma; myh++) {
                        sigma_off += nsopi_[myh];
                    }

                    int h_lambda_sigma = h_lambda ^ h_sigma;
                    if ( h_lambda_sigma != hia ) continue;

                    for (int lambda = 0; lambda < nsopi_[h_lambda]; lambda++) {
                        for (int sigma = 0; sigma < nsopi_[h_sigma]; sigma++) {

                            double dum_iajb = 0.0;
                            for (int mu = 0; mu < nsopi_[hi]; mu++) {
                                for (int nu = 0; nu < nsopi_[ha]; nu++) {
                                    dum_iajb += jelly_->ERI(mu+i_off,nu+a_off,lambda+lambda_off,sigma+sigma_off)
                                              * ci[mu    ][i            ]
                                              * ca[nu    ][a+doccpi_[ha]];

                                }
                            }
                            tmp_p[ia][(lambda+lambda_off)*nso_ + sigma+sigma_off] = dum_iajb;

                        }
                    }
                }
            }
        }
        // (ia|jb) second half transformation
        #pragma omp parallel for schedule(dynamic)
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** ci = Ca_->pointer(hi);
            double ** ca = Ca_->pointer(ha);

            for (int jb = 0; jb < cis_transition_list_.size(); jb++) {

                int j   = cis_transition_list_[jb].i;
                int b   = cis_transition_list_[jb].a;
                int hj  = cis_transition_list_[jb].hi;
                int hb  = cis_transition_list_[jb].ha;

                int j_off = 0;
                for (int myh = 0; myh < hj; myh++) {
                    j_off += nsopi_[myh];
                }

                int b_off = 0;
                for (int myh = 0; myh < hb; myh++) {
                    b_off += nsopi_[myh];
                }

                double ** cj = Ca_->pointer(hj);
                double ** cb = Ca_->pointer(hb);

                double dum_iajb = 0.0;
                for (int lambda = 0; lambda < nsopi_[hj]; lambda++) {
                    for (int sigma = 0; sigma < nsopi_[hb]; sigma++) {

                        dum_iajb += tmp_p[ia][(lambda+j_off)*nso_ + sigma+b_off]
                                  * cj[lambda][j            ]
                                  * cb[sigma ][b+doccpi_[hb]];
                    }
                }
                ham_p[ia][jb] = 2.0 * dum_iajb / Lfac_;
            }
        }
        tmp->zero();
        // (ij|ab) first half-transformation
        #pragma omp parallel for schedule(dynamic)
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;
            int hia = hi ^ ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** ci = Ca_->pointer(hi);
            double ** ca = Ca_->pointer(ha);

            for (int h_lambda = 0; h_lambda < nirrep_; h_lambda++) {
                int lambda_off = 0;
                for (int myh = 0; myh < h_lambda; myh++) {
                    lambda_off += nsopi_[myh];
                }
                for (int h_sigma = 0; h_sigma < nirrep_; h_sigma++) {
                    int sigma_off = 0;
                    for (int myh = 0; myh < h_sigma; myh++) {
                        sigma_off += nsopi_[myh];
                    }

                    int h_lambda_sigma = h_lambda ^ h_sigma;
                    if ( h_lambda_sigma != hia ) continue;

                    for (int lambda = 0; lambda < nsopi_[h_lambda]; lambda++) {
                        for (int sigma = 0; sigma < nsopi_[h_sigma]; sigma++) {

                            double dum_ijab = 0.0;
                            for (int mu = 0; mu < nsopi_[hi]; mu++) {
                                for (int nu = 0; nu < nsopi_[ha]; nu++) {
                                    dum_ijab += jelly_->ERI(mu+i_off,lambda+lambda_off,nu+a_off,sigma+sigma_off)
                                              * ci[mu    ][i            ]
                                              * ca[nu    ][a+doccpi_[ha]];

                                }
                            }
                            tmp_p[ia][(lambda+lambda_off)*nso_ + sigma+sigma_off] = dum_ijab;

                        }
                    }
                }
            }
        }
        // (ij|ab) second half transformation
        #pragma omp parallel for schedule(dynamic)
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** ci = Ca_->pointer(hi);
            double ** ca = Ca_->pointer(ha);

            for (int jb = 0; jb < cis_transition_list_.size(); jb++) {

                int j   = cis_transition_list_[jb].i;
                int b   = cis_transition_list_[jb].a;
                int hj  = cis_transition_list_[jb].hi;
                int hb  = cis_transition_list_[jb].ha;

                int j_off = 0;
                for (int myh = 0; myh < hj; myh++) {
                    j_off += nsopi_[myh];
                }

                int b_off = 0;
                for (int myh = 0; myh < hb; myh++) {
                    b_off += nsopi_[myh];
                }

                double ** cj = Ca_->pointer(hj);
                double ** cb = Ca_->pointer(hb);

                double dum_ijab = 0.0;
                for (int lambda = 0; lambda < nsopi_[hj]; lambda++) {
                    for (int sigma = 0; sigma < nsopi_[hb]; sigma++) {

                        dum_ijab += tmp_p[ia][(lambda+j_off)*nso_ + sigma+b_off]
                                  * cj[lambda][j            ]
                                  * cb[sigma ][b+doccpi_[hb]];
                    }
                }
                ham_p[ia][jb] -= dum_ijab / Lfac_;
            }
        }

        // remaining diagonal parts of Hamiltonian
        #pragma omp parallel for schedule(dynamic)
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** fa = F->pointer(ha);
            double ** fi = F->pointer(hi);
            ham_p[ia][ia] += fa[a+doccpi_[ha]][a+doccpi_[ha]] - fi[i][i];
        }
    }
    outfile->Printf("done.\n");

    outfile->Printf("    diagonalize CIS Hamiltonian..."); fflush(stdout);
    std::shared_ptr<Matrix> cis_eigvec (new Matrix(cis_ham));
    std::shared_ptr<Vector> cis_eigval = std::make_shared<Vector>(ovpi);
    cis_ham->diagonalize(cis_eigvec,cis_eigval,ascending);
    outfile->Printf("done.\n");

    outfile->Printf("\n");
    outfile->Printf("    ==> CIS excitaiton energies <=="); fflush(stdout);
    outfile->Printf("\n");

    for (int h = 0; h < nirrep_; h++) {

        int num_roots = 10;
        if ( options_["ROOTS_PER_IRREP"].has_changed() ) {
            num_roots = options_["ROOTS_PER_IRREP"][h].to_integer();
        }
        if ( num_roots == 0 ) continue;

        outfile->Printf("\n");
        outfile->Printf("    Irrep: %10s\n",jelly_->labels[h].c_str());
        outfile->Printf("\n");
        outfile->Printf("    state");
        outfile->Printf("                w(Eh)");
        outfile->Printf("                w(eV)");
        outfile->Printf("                    f");
        outfile->Printf("\n");
        double * eigval_p = cis_eigval->pointer(h);
        for (int I = 0; I < ovpi[h]; I++) {

            // evaluate oscillator strengths
            double tdp_x = sqrt(2.0) * C_DDOT(cis_transition_list_.size(),&(cis_eigvec->pointer(h)[0][I]),cis_transition_list_.size(),dipole_x->pointer(h),1);
            double tdp_y = sqrt(2.0) * C_DDOT(cis_transition_list_.size(),&(cis_eigvec->pointer(h)[0][I]),cis_transition_list_.size(),dipole_y->pointer(h),1);
            double tdp_z = sqrt(2.0) * C_DDOT(cis_transition_list_.size(),&(cis_eigvec->pointer(h)[0][I]),cis_transition_list_.size(),dipole_z->pointer(h),1);
            double f = 2.0 / 3.0 * eigval_p[I] * ( tdp_x * tdp_x + tdp_y * tdp_y + tdp_z * tdp_z );

            outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",I,eigval_p[I],eigval_p[I] * pc_hartree2ev, f);
        }
        outfile->Printf("\n");

    }

    // check sum
    double check = 0.0;
    for (int h = 0; h < nirrep_; h++) {
        check += C_DNRM2(ovpi[h],cis_eigval->pointer(h),1);
    }
    outfile->Printf("\n");
    outfile->Printf("    ||eps|| = %20.12lf\n",check);
    outfile->Printf("\n");

    //free(virpi);
    //free(ovpi);
}
void Jellium_SCFSolver::CIS_slow() {

    outfile->Printf("\n");
    outfile->Printf( "        ***************************************************\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    Finite Jellium CIS                           *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        ***************************************************\n");

    outfile->Printf("\n");
    int o = nelectron_ / 2;
    int v = nso_ - o;

    Dimension virpi(nirrep_);
    Dimension ovpi(nirrep_);
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
    std::shared_ptr<Matrix> F (new Matrix(Fa_));
    F->transform(Ca_);

    // dipole integrals
    std::shared_ptr<Vector> dipole_x = std::make_shared<Vector>(ovpi);
    std::shared_ptr<Vector> dipole_y = std::make_shared<Vector>(ovpi);
    std::shared_ptr<Vector> dipole_z = std::make_shared<Vector>(ovpi);

    /// transform ERIs to mo basis

    outfile->Printf("    transform (ia|jb), (ij|ab)...."); fflush(stdout);
    std::shared_ptr<Matrix> cis_ham (new Matrix(nirrep_,ovpi,ovpi));
    for (int h = 0; h < nirrep_; h++) {

        int num_roots = 10;
        if ( options_["ROOTS_PER_IRREP"].has_changed() ) {
            num_roots = options_["ROOTS_PER_IRREP"][h].to_integer();
        }
        if ( num_roots == 0 ) continue;

        if ( ovpi[h] == 0 ) continue;

        // list of CIS transitions in this irrep
        cis_transition_list_.clear();
        for (int hi = 0; hi < nirrep_; hi++) {
            for (int ha = 0; ha < nirrep_; ha++) {
                int hia = hi ^ ha;
                if ( hia != h ) continue;
                for (int i = 0; i < doccpi_[hi]; i++) {
                    for (int a = 0; a < virpi[ha]; a++) {

                        cis_transition me;
                        me.i   = i;
                        me.a   = a;
                        me.hi  = hi;
                        me.ha  = ha;

                        cis_transition_list_.push_back(me);
                    }
                }
            }
        }

        // dipole integrals
        double * dipole_x_p = dipole_x->pointer(h);
        double * dipole_y_p = dipole_y->pointer(h);
        double * dipole_z_p = dipole_z->pointer(h);
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** ci = Ca_->pointer(hi);
            double ** ca = Ca_->pointer(ha);

            double dum_x = 0.0;
            double dum_y = 0.0;
            double dum_z = 0.0;
            for (int mu = 0; mu < nsopi_[hi]; mu++) {
                for (int nu = 0; nu < nsopi_[ha]; nu++) {

                    double dip_x = jelly_->dipole_x(mu+i_off,nu+a_off,boxlength_);
                    double dip_y = jelly_->dipole_y(mu+i_off,nu+a_off,boxlength_);
                    double dip_z = jelly_->dipole_z(mu+i_off,nu+a_off,boxlength_);
                    dum_x += dip_x
                           * ci[mu][i            ]
                           * ca[nu][a+doccpi_[ha]];
                    dum_y += dip_y
                           * ci[mu][i            ]
                           * ca[nu][a+doccpi_[ha]];
                    dum_z += dip_z
                           * ci[mu][i            ]
                           * ca[nu][a+doccpi_[ha]];
                }
                dipole_x_p[ia] = dum_x;
                dipole_y_p[ia] = dum_y;
                dipole_z_p[ia] = dum_z;
            }
        }

        double ** ham_p = cis_ham->pointer(h);

        #pragma omp parallel for schedule(dynamic)
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** ci = Ca_->pointer(hi);
            double ** ca = Ca_->pointer(ha);

            for (int jb = 0; jb < cis_transition_list_.size(); jb++) {

                int j   = cis_transition_list_[jb].i;
                int b   = cis_transition_list_[jb].a;
                int hj  = cis_transition_list_[jb].hi;
                int hb  = cis_transition_list_[jb].ha;

                int j_off = 0;
                for (int myh = 0; myh < hj; myh++) {
                    j_off += nsopi_[myh];
                }

                int b_off = 0;
                for (int myh = 0; myh < hb; myh++) {
                    b_off += nsopi_[myh];
                }

                double ** cj = Ca_->pointer(hj);
                double ** cb = Ca_->pointer(hb);

                // (ia|jb) = (mu nu | lambda sigma) c(i,mu) c(a,nu) c(j,lambda) c(b,sigma)
                double dum_iajb = 0.0;
                for (int mu = 0; mu < nsopi_[hi]; mu++) {
                    for (int nu = 0; nu < nsopi_[ha]; nu++) {
                        for (int lambda = 0; lambda < nsopi_[hj]; lambda++) {
                            for (int sigma = 0; sigma < nsopi_[hb]; sigma++) {
                                dum_iajb += jelly_->ERI(mu+i_off,nu+a_off,lambda+j_off,sigma+b_off)
                                          * ci[mu    ][i            ]
                                          * ca[nu    ][a+doccpi_[ha]]
                                          * cj[lambda][j            ]
                                          * cb[sigma ][b+doccpi_[hb]];
                            }
                        }
                    }
                }

                // (ij|ab) = (mu nu | lambda sigma) c(i,mu) c(j,nu) c(a,lambda) c(b,sigma)
                double dum_ijab = 0.0;
                for (int mu = 0; mu < nsopi_[hi]; mu++) {
                    for (int nu = 0; nu < nsopi_[hj]; nu++) {
                        for (int lambda = 0; lambda < nsopi_[ha]; lambda++) {
                            for (int sigma = 0; sigma < nsopi_[hb]; sigma++) {
                                dum_ijab += jelly_->ERI(mu+i_off,nu+j_off,lambda+a_off,sigma+b_off)
                                          * ci[mu    ][i            ]
                                          * cj[nu    ][j            ]
                                          * ca[lambda][a+doccpi_[ha]]
                                          * cb[sigma ][b+doccpi_[hb]];
                            }
                        }
                    }
                }

                double dum = (2.0 * dum_iajb - dum_ijab)/Lfac_;
                if ( hi == hj && i == j ) {
                    double ** fp = F->pointer(ha);
                    dum += fp[a+doccpi_[ha]][b+doccpi_[ha]];
                }
                if ( ha == hb && a == b ) {
                    double ** fp = F->pointer(hi);
                    dum -= fp[i][j];
                }
                ham_p[ia][jb] = dum;
            }
        }
    }
    outfile->Printf("done.\n");

    outfile->Printf("    diagonalize CIS Hamiltonian..."); fflush(stdout);
    std::shared_ptr<Matrix> cis_eigvec (new Matrix(cis_ham));
    std::shared_ptr<Vector> cis_eigval = std::make_shared<Vector>(ovpi);
    cis_ham->diagonalize(cis_eigvec,cis_eigval,ascending);
    outfile->Printf("done.\n");

    outfile->Printf("\n");
    outfile->Printf("    ==> CIS excitaiton energies <=="); fflush(stdout);
    outfile->Printf("\n");

    for (int h = 0; h < nirrep_; h++) {

        int num_roots = 10;
        if ( options_["ROOTS_PER_IRREP"].has_changed() ) {
            num_roots = options_["ROOTS_PER_IRREP"][h].to_integer();
        }
        if ( num_roots == 0 ) continue;

        outfile->Printf("\n");
        outfile->Printf("    Irrep: %10s\n",jelly_->labels[h].c_str());
        outfile->Printf("\n");
        outfile->Printf("    state");
        outfile->Printf("                w(Eh)");
        outfile->Printf("                w(eV)");
        outfile->Printf("                    f");
        outfile->Printf("\n");
        double * eigval_p = cis_eigval->pointer(h);
        for (int I = 0; I < ovpi[h]; I++) {

            // evaluate oscillator strengths
            double tdp_x = sqrt(2.0) * C_DDOT(cis_transition_list_.size(),&(cis_eigvec->pointer(h)[0][I]),cis_transition_list_.size(),dipole_x->pointer(h),1);
            double tdp_y = sqrt(2.0) * C_DDOT(cis_transition_list_.size(),&(cis_eigvec->pointer(h)[0][I]),cis_transition_list_.size(),dipole_y->pointer(h),1);
            double tdp_z = sqrt(2.0) * C_DDOT(cis_transition_list_.size(),&(cis_eigvec->pointer(h)[0][I]),cis_transition_list_.size(),dipole_z->pointer(h),1);
            double f = 2.0 / 3.0 * eigval_p[I] * ( tdp_x * tdp_x + tdp_y * tdp_y + tdp_z * tdp_z );

            outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",I,eigval_p[I],eigval_p[I] * pc_hartree2ev, f);
        }
        outfile->Printf("\n");

    }

    // check sum
    double check = 0.0;
    for (int h = 0; h < nirrep_; h++) {
        check += C_DNRM2(ovpi[h],cis_eigval->pointer(h),1);
    }
    outfile->Printf("\n");
    outfile->Printf("    ||eps|| = %20.12lf\n",check);
    outfile->Printf("\n");

    //free(virpi);
    //free(ovpi);
}
void Jellium_SCFSolver::CIS_direct() {

    outfile->Printf("\n");
    outfile->Printf( "        ***************************************************\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *    Finite Jellium CIS                           *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        *                                                 *\n");
    outfile->Printf( "        ***************************************************\n");

    outfile->Printf("\n");
    int o = nelectron_ / 2;
    int v = nso_ - o;

    Dimension virpi(nirrep_);
    Dimension ovpi(nirrep_);
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

    // transform fock matrix to mo basis
    std::shared_ptr<Matrix> F (new Matrix(Fa_));
    F->transform(Ca_);

    // dipole integrals
    std::shared_ptr<Vector> dipole_x = std::make_shared<Vector>(ovpi);
    std::shared_ptr<Vector> dipole_y = std::make_shared<Vector>(ovpi);
    std::shared_ptr<Vector> dipole_z = std::make_shared<Vector>(ovpi);

    // diagonalize blocks by irrep
    double check = 0.0; 
    for (int h = 0; h < nirrep_; h++) {

        int num_roots = 10;
        if ( options_["ROOTS_PER_IRREP"].has_changed() ) {
            num_roots = options_["ROOTS_PER_IRREP"][h].to_integer();
        }
        if ( num_roots > ovpi[h] ) {
            num_roots = ovpi[h];
        }

        if ( num_roots == 0 ) continue;

        outfile->Printf("    Irrep:                     %10s\n",jelly_->labels[h].c_str());
        outfile->Printf("    Number of ia transitions:  %10i\n",ovpi[h]);
        outfile->Printf("    Number of roots requested: %10i\n",num_roots);
        outfile->Printf("\n");

        if ( ovpi[h] == 0 ) {
            outfile->Printf("\n");
            continue;
        }

        // list of CIS transitions in this irrep
        cis_transition_list_.clear();
        for (int hi = 0; hi < nirrep_; hi++) {
            for (int ha = 0; ha < nirrep_; ha++) {
                int hia = hi ^ ha;
                if ( hia != h ) continue;
                for (int i = 0; i < doccpi_[hi]; i++) {
                    for (int a = 0; a < virpi[ha]; a++) {

                        cis_transition me;
                        me.i   = i;
                        me.a   = a;
                        me.hi  = hi;
                        me.ha  = ha;

                        cis_transition_list_.push_back(me);
                    }
                }
            }
        }

        // dipole integrals
        outfile->Printf("    transform mu(i,a)................"); fflush(stdout);
        double * dipole_x_p = dipole_x->pointer(h);
        double * dipole_y_p = dipole_y->pointer(h);
        double * dipole_z_p = dipole_z->pointer(h);
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** ci = Ca_->pointer(hi);
            double ** ca = Ca_->pointer(ha);

            double dum_x = 0.0;
            double dum_y = 0.0;
            double dum_z = 0.0;
            for (int mu = 0; mu < nsopi_[hi]; mu++) {
                for (int nu = 0; nu < nsopi_[ha]; nu++) {

                    double dip_x = jelly_->dipole_x(mu+i_off,nu+a_off,boxlength_);
                    double dip_y = jelly_->dipole_y(mu+i_off,nu+a_off,boxlength_);
                    double dip_z = jelly_->dipole_z(mu+i_off,nu+a_off,boxlength_);
                    dum_x += dip_x
                           * ci[i            ][mu]
                           * ca[a+doccpi_[ha]][nu];
                    dum_y += dip_y
                           * ci[i            ][mu]
                           * ca[a+doccpi_[ha]][nu];
                    dum_z += dip_z
                           * ci[i            ][mu]
                           * ca[a+doccpi_[ha]][nu];
                }
                dipole_x_p[ia] = dum_x;
                dipole_y_p[ia] = dum_y;
                dipole_z_p[ia] = dum_z;
            }
        }
        outfile->Printf("done\n");

        // construct diagonal elements of CIS hamiltonian
        std::shared_ptr<Vector> cis_diagonal_ham (new Vector(cis_transition_list_.size()));
        double * ham_p = cis_diagonal_ham->pointer();
        // transform (ia|ia), (ii,aa)
        outfile->Printf("    transform (ia|ia), (ii|aa)......."); fflush(stdout);
        #pragma omp parallel for schedule(dynamic)
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** ci = Ca_->pointer(hi);
            double ** ca = Ca_->pointer(ha);

            // (ia|ia) = (mu nu | lambda sigma) c(i,mu) c(a,nu) c(i,lambda) c(a,sigma)
            double dum_iaia = 0.0;
            for (int mu = 0; mu < nsopi_[hi]; mu++) {
                for (int nu = 0; nu < nsopi_[ha]; nu++) {

                    for (int lambda = 0; lambda < nsopi_[hi]; lambda++) {
                        for (int sigma = 0; sigma < nsopi_[ha]; sigma++) {
                            double eri = jelly_->ERI(mu+i_off,nu+a_off,lambda+i_off,sigma+a_off);
                            dum_iaia += eri 
                                      * ci[i            ][mu]
                                      * ca[a+doccpi_[ha]][nu]
                                      * ci[i            ][lambda]
                                      * ca[a+doccpi_[ha]][sigma];
                        }
                    }
                }
            }

            // (ii|aa) = (mu nu | lambda sigma) c(i,mu) c(i,nu) c(a,lambda) c(a,sigma)
            double dum_iiaa = 0.0;
            for (int mu = 0; mu < nsopi_[hi]; mu++) {
                for (int nu = 0; nu < nsopi_[hi]; nu++) {
                    for (int lambda = 0; lambda < nsopi_[ha]; lambda++) {
                        for (int sigma = 0; sigma < nsopi_[ha]; sigma++) {
                            double eri = jelly_->ERI(mu+i_off,nu+i_off,lambda+a_off,sigma+a_off);
                            dum_iiaa += eri
                                      * ci[i            ][mu]
                                      * ci[i            ][nu]
                                      * ca[a+doccpi_[ha]][lambda]
                                      * ca[a+doccpi_[ha]][sigma];
                        }
                    }
                }
            }

            double dum = (2.0 * dum_iaia - dum_iiaa)/Lfac_;
            dum += F->pointer(ha)[a+doccpi_[ha]][a+doccpi_[ha]];
            dum -= F->pointer(hi)[i][i];
            
            ham_p[ia] = dum;

        }
        outfile->Printf("done\n");

        outfile->Printf("    diagonalize CIS Hamiltonian......"); fflush(stdout);
        std::shared_ptr<DavidsonSolver> david (new DavidsonSolver());

        davidson_irrep_ = h;

        size_t ci_iter = 0;
        int print      = 1;

        std::shared_ptr<Matrix> eigvec (new Matrix(num_roots,ovpi[h]));
        std::shared_ptr<Vector> eigval (new Vector(num_roots));

        // type of guess for davidson. the smart guess will do N^4 work for each hamiltonian element
        HamiltonianElementFunction my_hamiltonian_element;
        if ( options_.get_bool("JELLIUM_CIS_SMART_GUESS") ) {
            my_hamiltonian_element = hamiltonian_element;
        }else {
            my_hamiltonian_element = NULL;
        }

        david->solve(cis_diagonal_ham->pointer(),
            ovpi[h],
            num_roots,
            eigval->pointer(),
            eigvec->pointer(),
            options_.get_double("R_CONVERGENCE"),
            print,
            evaluate_sigma,
            my_hamiltonian_element,
            ci_iter,
            (void*)this,
            options_.get_int("DAVIDSON_MAXDIM") * num_roots);
        outfile->Printf("done\n");

        outfile->Printf("    CIS excitation energies:\n");
        outfile->Printf("\n");
        outfile->Printf("    state");
        outfile->Printf("                w(Eh)");
        outfile->Printf("                w(eV)");
        outfile->Printf("                    f");
        outfile->Printf("\n");
        double * eigval_p = eigval->pointer();
        for (int I = 0; I < num_roots; I++) {

            // evaluate oscillator strengths
            double tdp_x = sqrt(2.0) * C_DDOT(cis_transition_list_.size(),eigvec->pointer()[I],1,dipole_x->pointer(h),1);
            double tdp_y = sqrt(2.0) * C_DDOT(cis_transition_list_.size(),eigvec->pointer()[I],1,dipole_y->pointer(h),1);
            double tdp_z = sqrt(2.0) * C_DDOT(cis_transition_list_.size(),eigvec->pointer()[I],1,dipole_z->pointer(h),1);
            double f = 2.0 / 3.0 * eigval_p[I] * ( tdp_x * tdp_x + tdp_y * tdp_y + tdp_z * tdp_z );

            outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",I,eigval_p[I],eigval_p[I] * pc_hartree2ev, f);
        }
        outfile->Printf("\n");

        check += C_DNRM2(ovpi[h],eigval->pointer(),1);
    }

    // check sum
    outfile->Printf("\n");
    outfile->Printf("    ||eps|| = %20.12lf\n",check);
    outfile->Printf("\n");

    //free(virpi);
    //free(ovpi);
}

void Jellium_SCFSolver::CIS_evaluate_sigma(size_t N, size_t maxdim, double ** bmat, double ** sigma) {

    for (size_t ia = 0; ia < N; ia++) {
        for (size_t k = 0; k < maxdim; k++) {
            sigma[ia][k] = 0.0;
        }
    }

    /// transform fock matrix to mo basis
    std::shared_ptr<Matrix> F (new Matrix(Fa_));
    F->transform(Ca_);

    //int * jb = (int*)malloc(nirrep_*sizeof(int));

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

    //  sigma[ ia ][ k ] = Ham[ ia ][ jb] *  bmat[ k ][ jb ], k = a given sigma vector, ia and jb = CIS configurations
    for (int k = 0; k < maxdim; k++) {

        // build dressed density matrix D(mu,nu) = Ca(j,mu) Ca(b,nu) b[k][jb], which i guess has no symmetry

        std::shared_ptr<Matrix> D (new Matrix(nso_,nso_));
        double ** dp = D->pointer();

        // loop over all CIS transitions, jb
        for (int jb = 0; jb < cis_transition_list_.size(); jb++) {

            int j   = cis_transition_list_[jb].i;
            int b   = cis_transition_list_[jb].a;
            int hj  = cis_transition_list_[jb].hi;
            int hb  = cis_transition_list_[jb].ha;

            int mu_off = 0;
            for (int myh = 0; myh < hj; myh++) {
                mu_off += nsopi_[myh];
            }

            int nu_off = 0;
            for (int myh = 0; myh < hb; myh++) {
                nu_off += nsopi_[myh];
            }

            double ** cj = Ca_->pointer(hj);
            double ** cb = Ca_->pointer(hb);

            for (int mu = 0; mu < nsopi_[hj]; mu++) {
                for (int nu = 0; nu < nsopi_[hb]; nu++) {
                    dp[mu + mu_off][nu + nu_off] += cj[mu][j] * cb[nu][b+doccpi_[hb]] * bmat[k][jb];
                }
            }
        }

        // build dressed coulomb matrix J(mu,nu) = D(lambda,sigma) (mu nu| lambda sigma), which i guess has no symmetry
        // build dressed exchange matrix K(mu,lambda) = D(nu,sigma) (mu nu| lambda sigma), which i guess has no symmetry

        std::shared_ptr<Matrix> Jmat (new Matrix(nso_,nso_));
        std::shared_ptr<Matrix> Kmat (new Matrix(nso_,nso_));

        double ** jp = Jmat->pointer();
        double ** kp = Kmat->pointer();

        for (int h_mu = 0; h_mu < nirrep_; h_mu++) {

            int mu_off = 0;
            for (int myh = 0; myh < h_mu; myh++) {
                mu_off += nsopi_[myh];
            }

            for (int h_nu = 0; h_nu < nirrep_; h_nu++) {

                int nu_off = 0;
                for (int myh = 0; myh < h_nu; myh++) {
                    nu_off += nsopi_[myh];
                }

                int h_mu_nu = h_mu ^ h_nu;

                #pragma omp parallel for schedule(dynamic)
                for (int mu = 0; mu < nsopi_[h_mu]; mu++) {
                    for (int nu = 0; nu < nsopi_[h_nu]; nu++) {

                        double myJ = 0.0;

                        for (int h_lambda = 0; h_lambda < nirrep_; h_lambda++) {

                            int lambda_off = 0;
                            for (int myh = 0; myh < h_lambda; myh++) {
                                lambda_off += nsopi_[myh];
                            }

                            for (int h_sigma = 0; h_sigma < nirrep_; h_sigma++) {

                                int sigma_off = 0;
                                for (int myh = 0; myh < h_sigma; myh++) {
                                    sigma_off += nsopi_[myh];
                                }

                                int h_lambda_sigma = h_lambda ^ h_sigma;
                                if ( h_mu_nu != h_lambda_sigma ) continue;

                                for (int lambda = 0; lambda < nsopi_[h_lambda]; lambda++) {
                                    for (int sigma = 0; sigma < nsopi_[h_sigma]; sigma++) {
                                        double dJ  = dp[lambda+lambda_off][sigma+sigma_off];
                                        double dK  = dp[nu+nu_off][sigma+sigma_off];
                                        double eri = jelly_->ERI(mu+mu_off,nu+nu_off,lambda+lambda_off,sigma+sigma_off);
                                        myJ += dJ * eri;
                                        kp[mu+mu_off][lambda+lambda_off] += dK * eri;
                                    }
                                }
                            }
                        }
                        jp[mu+mu_off][nu+nu_off] = myJ;
                    }
                }
            }
        }
        Jmat->scale(1.0/Lfac_);
        Kmat->scale(1.0/Lfac_);

        // loop over all CIS transitions, ia
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int mu_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                mu_off += nsopi_[myh];
            }

            int nu_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                nu_off += nsopi_[myh];
            }

            double ** ci = Ca_->pointer(hi);
            double ** ca = Ca_->pointer(ha);

            // fock matrix contribution on diagonal
            double dum = 0.0;
            dum += F->pointer(ha)[a+doccpi_[ha]][a+doccpi_[ha]];
            dum -= F->pointer(hi)[i][i];

            sigma[ia][k] = dum * bmat[k][ia];

            // 2 (ia|jb) - (ij|ab) contribution:
            // [ 2 J(mu,nu) - K(mu,nu) ] c(i,mu) c(a,nu)
            
            dum = 0.0;
            for (int mu = 0; mu < nsopi_[hi]; mu++) {
                for (int nu = 0; nu < nsopi_[ha]; nu++) {
                    dum += (2.0 * jp[mu+mu_off][nu+nu_off] - kp[mu+mu_off][nu+nu_off]) 
                         * ci[mu][i            ] 
                         * ca[nu][a+doccpi_[ha]];
                }
            }

            sigma[ia][k] += dum;
        }
    }

    //free(jb);
}

double Jellium_SCFSolver::evaluate_hamiltonian_element(size_t my_ia, size_t my_jb) {

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
    std::shared_ptr<Matrix> F (new Matrix(Fa_));
    F->transform(Ca_);

    /// transform ERIs to mo basis

    //outfile->Printf("    transform (ia|jb), (ij|ab)...."); fflush(stdout);
    //std::shared_ptr<Matrix> cis_ham (new Matrix(nirrep_,ovpi,ovpi));
    double return_value;
    for (int h = 0; h < nirrep_; h++) {

        // clunky, but works
        if ( h != davidson_irrep_ ) continue;

        if ( ovpi[h] == 0 ) continue;

        // list of CIS transitions in this irrep
        cis_transition_list_.clear();
        for (int hi = 0; hi < nirrep_; hi++) {
            for (int ha = 0; ha < nirrep_; ha++) {
                int hia = hi ^ ha;
                if ( hia != h ) continue;
                for (int i = 0; i < doccpi_[hi]; i++) {
                    for (int a = 0; a < virpi[ha]; a++) {

                        cis_transition me;
                        me.i   = i;
                        me.a   = a;
                        me.hi  = hi;
                        me.ha  = ha;

                        cis_transition_list_.push_back(me);
                    }
                }
            }
        }

        //double ** ham_p = cis_ham->pointer(h);

        #pragma omp parallel for schedule(dynamic)
        for (int ia = 0; ia < cis_transition_list_.size(); ia++) {

            if ( ia != my_ia ) continue;

            int i   = cis_transition_list_[ia].i;
            int a   = cis_transition_list_[ia].a;
            int hi  = cis_transition_list_[ia].hi;
            int ha  = cis_transition_list_[ia].ha;

            int i_off = 0;
            for (int myh = 0; myh < hi; myh++) {
                i_off += nsopi_[myh];
            }

            int a_off = 0;
            for (int myh = 0; myh < ha; myh++) {
                a_off += nsopi_[myh];
            }

            double ** ci = Ca_->pointer(hi);
            double ** ca = Ca_->pointer(ha);

            for (int jb = 0; jb < cis_transition_list_.size(); jb++) {

                if ( jb != my_jb ) continue;

                int j   = cis_transition_list_[jb].i;
                int b   = cis_transition_list_[jb].a;
                int hj  = cis_transition_list_[jb].hi;
                int hb  = cis_transition_list_[jb].ha;

                int j_off = 0;
                for (int myh = 0; myh < hj; myh++) {
                    j_off += nsopi_[myh];
                }

                int b_off = 0;
                for (int myh = 0; myh < hb; myh++) {
                    b_off += nsopi_[myh];
                }

                double ** cj = Ca_->pointer(hj);
                double ** cb = Ca_->pointer(hb);

                // (ia|jb) = (mu nu | lambda sigma) c(i,mu) c(a,nu) c(j,lambda) c(b,sigma)
                double dum_iajb = 0.0;
                for (int mu = 0; mu < nsopi_[hi]; mu++) {
                    for (int nu = 0; nu < nsopi_[ha]; nu++) {
                        for (int lambda = 0; lambda < nsopi_[hj]; lambda++) {
                            for (int sigma = 0; sigma < nsopi_[hb]; sigma++) {
                                dum_iajb += jelly_->ERI(mu+i_off,nu+a_off,lambda+j_off,sigma+b_off)
                                          * ci[i            ][mu]
                                          * ca[a+doccpi_[ha]][nu]
                                          * cj[j            ][lambda]
                                          * cb[b+doccpi_[hb]][sigma];
                            }
                        }
                    }
                }

                // (ij|ab) = (mu nu | lambda sigma) c(i,mu) c(j,nu) c(a,lambda) c(b,sigma)
                double dum_ijab = 0.0;
                for (int mu = 0; mu < nsopi_[hi]; mu++) {
                    for (int nu = 0; nu < nsopi_[hj]; nu++) {
                        for (int lambda = 0; lambda < nsopi_[ha]; lambda++) {
                            for (int sigma = 0; sigma < nsopi_[hb]; sigma++) {
                                dum_ijab += jelly_->ERI(mu+i_off,nu+j_off,lambda+a_off,sigma+b_off)
                                          * ci[i            ][mu]
                                          * cj[j            ][nu]
                                          * ca[a+doccpi_[ha]][lambda]
                                          * cb[b+doccpi_[hb]][sigma];
                            }
                        }
                    }
                }

                double dum = (2.0 * dum_iajb - dum_ijab)/Lfac_;
                if ( hi == hj && i == j ) {
                    double ** fp = F->pointer(ha);
                    dum += fp[a+doccpi_[ha]][b+doccpi_[ha]];
                }
                if ( ha == hb && a == b ) {
                    double ** fp = F->pointer(hi);
                    dum -= fp[i][j];
                }
                return_value =  dum;
                break;
                //ham_p[ia][jb] = dum;
            }
        }
    }
    //outfile->Printf("done.\n");
    free(virpi);
    free(ovpi);
    return return_value;

}


} // End namespaces



