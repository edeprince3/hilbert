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
}

Jellium_SCFSolver::~Jellium_SCFSolver(){
}

double Jellium_SCFSolver::compute_energy(){

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
    std::shared_ptr<JelliumIntegrals> jelly (new JelliumIntegrals(options_));

    // paramters
    int nso       = options_.get_int("N_BASIS_FUNCTIONS"); 
    int nelectron = options_.get_int("N_ELECTRONS");
    int * nsopi   = jelly->nsopi();
    int * doccpi  = jelly->doccpi();
    int nirrep    = jelly->nirrep();

    // number of electrons
    if ( nelectron % 2 != 0 ) {
        throw PsiException("jellium currently requires an even number of electrons.",__FILE__,__LINE__);
    }

    // factor for box size ... coded to give <rho> = 1
    double Lfac = pow(nelectron,1.0/3.0) / M_PI;
    double boxlength = M_PI / 2; 
    
    //since box is already pi a.u. long
    //double length_nm = options_.get_double("length");
    //Lfac = length_nm * 18.89725988 / M_PI;
    //boxlength = options_.get_double("LENGTH");
    
    // kinetic energy integrals
    std::shared_ptr<Matrix> T = jelly->Ke;
    T->scale(1.0/Lfac/Lfac);

    // potential energy integrals
    std::shared_ptr<Matrix> V = jelly->NucAttrac;
    V->scale(1.0/Lfac);

    // convergence parameters
    double e_convergence = options_.get_double("E_CONVERGENCE");
    double d_convergence = options_.get_double("D_CONVERGENCE");
    int    maxiter       = options_.get_int("MAXITER");
    bool   do_diis       = options_.get_bool("DIIS");

    // print some information about this computation
    outfile->Printf("\n");
    outfile->Printf("    ==> Jellium Hartree-Fock parameters <==\n");
    outfile->Printf("\n");
    outfile->Printf("    number of electrons:              %10i\n",nelectron);
    outfile->Printf("    number of basis functions:        %10i\n",nso);
    outfile->Printf("    maximum particle-in-a-box state:  %10i\n",jelly->get_nmax());
    outfile->Printf("    e_convergence                     %10.2le\n",e_convergence);
    outfile->Printf("    d_convergence                     %10.2le\n",d_convergence);
    outfile->Printf("    maxiter                           %10i\n",maxiter);
    outfile->Printf("    diis?                             %10s\n",do_diis ? "yes" : "no");
    outfile->Printf("\n");

    std::shared_ptr<Matrix> Ca (new Matrix(nirrep,nsopi,nsopi));

    std::shared_ptr<DIIS> diis (new DIIS(nso*nso));
    
    // build core hamiltonian
    std::shared_ptr<Matrix> h (new Matrix(T));
    V->scale(nelectron); 
    h->add(V);

    // fock matrix
    std::shared_ptr<Matrix> Fa = (std::shared_ptr<Matrix>)(new Matrix(h));

    // eigenvectors / eigenvalues of fock matrix
    std::shared_ptr<Vector> Feval (new Vector(nirrep,nsopi));

    // diagonalize core hamiltonian, get orbitals
    Fa->diagonalize(Ca,Feval);

    // build density matrix 
    std::shared_ptr<Matrix> Da (new Matrix(nirrep,nsopi,nsopi));
    for (int h = 0; h < nirrep; h++) {
        C_DGEMM('n','t',nsopi[h],nsopi[h],doccpi[h],1.0,
            &(Ca->pointer(h)[0][0]),nsopi[h],
            &(Ca->pointer(h)[0][0]),nsopi[h],0.0,
            &(Da->pointer(h)[0][0]),nsopi[h]);
    }

    double energy = Da->vector_dot(h) + Da->vector_dot(Fa);
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
   
    do {

        // containers for J and K 
        std::shared_ptr<Matrix> Ja (new Matrix(nirrep,nsopi,nsopi));
        std::shared_ptr<Matrix> Ka (new Matrix(nirrep,nsopi,nsopi));

        for (int hp = 0; hp < nirrep; hp++) {
            int offp = 0;
            double ** k_p = Ka->pointer(hp);
            double ** j_p = Ja->pointer(hp);
            for (int myh = 0; myh < hp; myh++) {
                offp += nsopi[myh];
            }

            #pragma omp parallel for schedule(dynamic)
            for (int p = 0; p < nsopi[hp]; p++) {
                int pp = p + offp;

                for (int q = p; q < nsopi[hp]; q++) {
                    int qq = q + offp;
                    double myJ = 0.0;
                    double myK = 0.0;
                    
                    for (int hr = 0; hr < nirrep; hr++) {
                        double ** d_p = Da->pointer(hr);
                        int offr = 0;
                        
                        for (int myh = 0; myh < hr; myh++) {
                            offr += nsopi[myh];
                        }

                        // exchange
                        for (int r = 0; r < nsopi[hr]; r++) {
                            int rr = r + offr;
                            for (int s = 0; s < nsopi[hr]; s++) {
                                int ss = s + offr;
                                myK += d_p[r][s] * jelly->ERI_int(pp,ss,rr,qq);
                            }
                        }

                        // coulomb
                        for (int r = 0; r < nsopi[hr]; r++) {
                            int rr = r + offr;
                            double dum = 0.0;
                            for (int s = r+1; s < nsopi[hr]; s++) {
                                int ss = s + offr;
                                dum += d_p[r][s] * jelly->ERI_int(pp,qq,rr,ss);
                            }
                            myJ += 2.0 * dum + d_p[r][r] * jelly->ERI_int(pp,qq,rr,rr);
                        }
                    }
                    j_p[p][q] = myJ;
                    j_p[q][p] = myJ;
                    k_p[p][q] = myK;
                    k_p[q][p] = myK;
                }
            }
        }
        if(iter>=2){jelly->iter = 2;}

        // fock matrix
        Fa->copy(Ja);
        Fa->scale(2.0);
        Fa->subtract(Ka);
        Fa->scale(1.0/Lfac);
        Fa->add(h);

        // evaluate current energy, E = D(H+F)
        double new_energy = 0.0;
        new_energy += (nelectron*nelectron/2.0) * jelly->selfval/Lfac;
        new_energy += Da->vector_dot(h);
        new_energy += Da->vector_dot(Fa);

        // orbital gradient
        std::shared_ptr<Matrix> orbital_gradient(new Matrix("[F,D]", nirrep, nsopi, nsopi));
        orbital_gradient->gemm(false, false, 1.0, Fa, Da, 0.0);

        std::shared_ptr<Matrix> tmp(orbital_gradient->transpose());
        orbital_gradient->subtract(tmp);
        tmp.reset();

        g_nrm = 0.0;
        for (int h = 0; h < nirrep; h++) {
            double dum = C_DNRM2(nsopi[h]*nsopi[h],orbital_gradient->pointer(h)[0],1);;
            g_nrm += dum * dum;
        }
        g_nrm = sqrt(g_nrm);

        // DIIS
        if (do_diis) {

            // diis wants data in vector
           
            std::shared_ptr<Vector> tmp_F(new Vector(nso*nso));
            std::shared_ptr<Vector> tmp_G(new Vector(nso*nso));

            int off = 0;
            for (int h = 0; h < nirrep; h++) {
                C_DCOPY(nsopi[h]*nsopi[h],Fa->pointer(h)[0],1,tmp_F->pointer() + off,1);
                C_DCOPY(nsopi[h]*nsopi[h],orbital_gradient->pointer(h)[0],1,tmp_G->pointer() + off,1);
                off += nsopi[h]*nsopi[h];
            }

            diis->WriteVector(tmp_F->pointer());
            diis->WriteErrorVector(tmp_G->pointer());
            diis->Extrapolate(tmp_F->pointer());

            off = 0;
            for (int h = 0; h < nirrep; h++) {
                C_DCOPY(nsopi[h]*nsopi[h],tmp_F->pointer()+off,1,Fa->pointer(h)[0],1);
                C_DCOPY(nsopi[h]*nsopi[h],tmp_G->pointer() + off,1,orbital_gradient->pointer(h)[0],1);
                off += nsopi[h]*nsopi[h];
            }

        }

        dele = fabs(new_energy - energy);

        // evaluate new orbitals
        std::shared_ptr<Vector> Feval (new Vector(nirrep, nsopi));
        Fa->diagonalize(Ca,Feval);

        // evaluate new density
        for (int h = 0; h < nirrep; h++) {
            C_DGEMM('n','t',nsopi[h],nsopi[h],doccpi[h],1.0,
                &(Ca->pointer(h)[0][0]),nsopi[h],
                &(Ca->pointer(h)[0][0]),nsopi[h],0.0,
                &(Da->pointer(h)[0][0]),nsopi[h]);
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

    //double fock_energy = D->vector_dot(K) / Lfac;
    //outfile->Printf("      Fock energy:             %20.12lf\n", fock_energy);

    outfile->Printf("    * Jellium SCF total energy: %20.12lf\n", energy);


    Process::environment.globals["CURRENT ENERGY"]    = energy;
    Process::environment.globals["JELLIUM SCF TOTAL ENERGY"] = energy;

    return energy;

}

} // End namespaces



