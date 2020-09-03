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

// jk object
#include <psi4/libfock/jk.h>

// diis solver
#include <misc/diis.h>

#include "rhf.h"

namespace hilbert{ 

PolaritonicRHF::PolaritonicRHF(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_):
    PolaritonicHF(reference_wavefunction,options_) {
    common_init();
}

PolaritonicRHF::~PolaritonicRHF() {
}

void PolaritonicRHF::common_init() {

    outfile->Printf("\n\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Polaritonic RHF                                  *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");

    // ensure scf_type df
    if ( options_.get_str("SCF_TYPE") != "DF" ) {
        throw PsiException("polaritonic rhf only works with scf_type df for now",__FILE__,__LINE__);
    }

    // ensure running in c1 symmetry
    if ( reference_wavefunction_->nirrep() > 1 ) {
        throw PsiException("polaritonic rhf only works with c1 symmetry for now.",__FILE__,__LINE__);
    }

    // ensure closed shell
    if ( nalpha_ != nbeta_ ) {
        throw PsiException("polaritonic rhf only works for closed shells",__FILE__,__LINE__);
    }

}

double PolaritonicRHF::compute_energy() {

    // grab the one-electron integrals from MintsHelper:
    std::shared_ptr<MintsHelper> mints (new MintsHelper(reference_wavefunction_));

    // one-electron kinetic energy integrals
    std::shared_ptr<Matrix> T = mints->so_kinetic();

    // one-electron potential energy integrals
    std::shared_ptr<Matrix> V = mints->so_potential();

    // build the core hamiltonian
    std::shared_ptr<Matrix> h = (std::shared_ptr<Matrix>)(new Matrix(T));
    h->add(V);

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    // get auxiliary basis:
    std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

    // total number of auxiliary basis functions
    int nQ = auxiliary->nbf();

    // JK object
    std::shared_ptr<DiskDFJK> jk = (std::shared_ptr<DiskDFJK>)(new DiskDFJK(primary,auxiliary));

    // memory for jk (say, 80% of what is available)
    jk->set_memory(0.8 * memory_);

    // integral cutoff
    jk->set_cutoff(options_.get_double("INTS_TOLERANCE"));

    // Do J/K, Not wK 
    jk->set_do_J(true);
    jk->set_do_K(true);
    jk->set_do_wK(false);

    jk->initialize();

    // grab some input options_
    double e_convergence = options_.get_double("E_CONVERGENCE");
    double d_convergence = options_.get_double("D_CONVERGENCE");
    int maxiter          = options_.get_int("MAXITER");

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso_);
    outfile->Printf("    No. auxiliary basis functions:  %5i\n",nQ);
    outfile->Printf("    No. electrons:                  %5i\n",nalpha_ + nbeta_);
    outfile->Printf("    e_convergence:             %10.3le\n",e_convergence);
    outfile->Printf("    d_convergence:             %10.3le\n",d_convergence);
    outfile->Printf("    maxiter:                        %5i\n",maxiter);
    outfile->Printf("\n");
    outfile->Printf("\n");

    // allocate memory for eigenvectors and eigenvalues of the overlap matrix
    std::shared_ptr<Matrix> Sevec ( new Matrix(nso_,nso_) );
    std::shared_ptr<Vector> Seval ( new Vector(nso_) );

    // build S^(-1/2) symmetric orthogonalization matrix
    S_->diagonalize(Sevec,Seval);

    std::shared_ptr<Matrix> Shalf = (std::shared_ptr<Matrix>)( new Matrix(nso_,nso_) );
    for (int mu = 0; mu < nso_; mu++) {
        Shalf->pointer()[mu][mu] = 1.0 / sqrt(Seval->pointer()[mu]);
    }

    // transform Seval back to nonorthogonal basis
    Shalf->back_transform(Sevec);

    // allocate memory for F' and its eigenvectors and eigenvalues
    std::shared_ptr<Matrix> Fevec ( new Matrix(nso_,nso_) );
    std::shared_ptr<Vector> Feval ( new Vector("Orbital Energies",nso_) );
    std::shared_ptr<Matrix> Fprime ( new Matrix(Fa_) );

/*
    // core guess ... shouldn't be necessary since we're starting from an existing reference

    // form F' = ST^(-1/2) F S^(-1/2), where F = h
    Fa_->copy(h);
    Fprime->transform(Shalf);

    // diagonalize F' to obtain C'
    Fprime->diagonalize(Fevec,Feval,ascending);

    // Find C = S^(-1/2)C'
    Ca_->gemm(false,false,1.0,Shalf,Fevec,0.0);

    // Construct density from C
    C_DGEMM('n','t',nso_,nso_,nalpha_,1.0,&(Ca_->pointer()[0][0]),nso_,&(Ca_->pointer()[0][0]),nso_,0.0,&(Fprime->pointer()[0][0]),nso_);
    Da_->copy(Fprime);
*/

    energy_  = enuc_;
    energy_ += Da_->vector_dot(h);
    energy_ += Da_->vector_dot(Fa_);

    // SCF iterations

    double e_last    = 0.0;
    double dele      = 0.0;
    double gnorm     = 0.0;

    outfile->Printf("\n");
    outfile->Printf("    Guess energy:  %20.12lf\n",energy_);
    outfile->Printf("\n");
    outfile->Printf("    ==>  Begin SCF Iterations <==\n");
    outfile->Printf("\n");
    outfile->Printf("    ");
    outfile->Printf(" Iter ");
    outfile->Printf("              energy ");
    outfile->Printf("                  dE ");
    outfile->Printf("         RMS |[F,P]| ");
    outfile->Printf("\n");

    std::shared_ptr<DIIS> diis (new DIIS(nso_*nso_));

    int iter = 0;
    do {

        e_last = energy_;

        // grab occupied orbitals (the first nalpha)
        std::shared_ptr<Matrix> myC (new Matrix(Ca_) );
        myC->zero();

        for (int mu = 0; mu < nso_; mu++) {
            for (int i = 0; i < nalpha_; i++) {
                myC->pointer()[mu][i] = Ca_->pointer()[mu][i];
            }
        }

        // push occupied orbitals onto JK object
        std::vector< std::shared_ptr<Matrix> >& C_left  = jk->C_left();
        C_left.clear();
        C_left.push_back(myC);

        // form J/K
        jk->compute();

        // form F = h + 2*J - K 
        Fa_->copy(jk->J()[0]);
        Fa_->scale(2.0);
        Fa_->subtract(jk->K()[0]);
        Fa_->add(h);

        if ( options_.get_int("N_PHOTON_STATES") > 1 ) {

            build_cavity_hamiltonian();

            std::shared_ptr<Matrix> Vx = (std::shared_ptr<Matrix>)(new Matrix(CavityDipolePotential_x_));
            std::shared_ptr<Matrix> Vy = (std::shared_ptr<Matrix>)(new Matrix(CavityDipolePotential_y_));
            std::shared_ptr<Matrix> Vz = (std::shared_ptr<Matrix>)(new Matrix(CavityDipolePotential_z_));

            Vx->scale(-CavityDipole_x_->pointer()[0][0]);
            Vy->scale(-CavityDipole_y_->pointer()[0][0]);
            Vz->scale(-CavityDipole_z_->pointer()[0][0]);

            Fa_->add(Vx);
            Fa_->add(Vy);
            Fa_->add(Vz);

        }

        // Construct density from C
        C_DGEMM('n','t',nso_,nso_,nalpha_,1.0,&(Ca_->pointer()[0][0]),nso_,&(Ca_->pointer()[0][0]),nso_,0.0,&(Da_->pointer()[0][0]),nso_);

        // evaluate the current energy, E = D(H+F) + Enuc
        energy_  = enuc_;
        energy_ += Da_->vector_dot(h);
        energy_ += Da_->vector_dot(Fa_);

        if (options_.get_int("N_PHOTON_STATES") > 1) {

            build_cavity_hamiltonian();

            std::shared_ptr<Matrix> Vx = (std::shared_ptr<Matrix>)(new Matrix(CavityDipolePotential_x_));
            std::shared_ptr<Matrix> Vy = (std::shared_ptr<Matrix>)(new Matrix(CavityDipolePotential_y_));
            std::shared_ptr<Matrix> Vz = (std::shared_ptr<Matrix>)(new Matrix(CavityDipolePotential_z_));

            Vx->scale(-CavityDipole_x_->pointer()[0][0]);
            Vy->scale(-CavityDipole_y_->pointer()[0][0]);
            Vz->scale(-CavityDipole_z_->pointer()[0][0]);

            energy_ += 2.0 * Da_->vector_dot(Vx);
            energy_ += 2.0 * Da_->vector_dot(Vy);
            energy_ += 2.0 * Da_->vector_dot(Vz);

        }

        // dele
        dele = energy_ - e_last;

        // form F' = ST^(-1/2) F S^(-1/2)
        Fprime->copy(Fa_);
        Fprime->transform(Shalf);

        // The error vector in DIIS for SCF is defined as 
        // the orbital gradient, in the orthonormal basis:
        // 
        // ST^{-1/2} [FDS - SDF] S^{-1/2}
        
        std::shared_ptr<Matrix> FDSmSDF(new Matrix("FDS-SDF", nso_, nso_));
        std::shared_ptr<Matrix> DS(new Matrix("DS", nso_, nso_));

        DS->gemm(false,false,1.0,Da_,S_,0.0);
        FDSmSDF->gemm(false,false,1.0,Fa_,DS,0.0);

        DS.reset();

        std::shared_ptr<Matrix> SDF(FDSmSDF->transpose());
        FDSmSDF->subtract(SDF);

        SDF.reset();

        std::shared_ptr<Matrix> ShalfGrad(new Matrix("ST^{-1/2}(FDS - SDF)", nso_, nso_));
        ShalfGrad->gemm(true,false,1.0,Shalf,FDSmSDF,0.0);
        FDSmSDF.reset();

        std::shared_ptr<Matrix> ShalfGradShalf(new Matrix("ST^{-1/2}(FDS - SDF)S^{-1/2}", nso_, nso_));
        ShalfGradShalf->gemm(false,false,1.0,ShalfGrad,Shalf,0.0);

        ShalfGrad.reset();

        // We will use the RMS of the orbital gradient 
        // to monitor convergence.
        gnorm = ShalfGradShalf->rms();

        // DIIS extrapolation
        diis->WriteVector(&(Fprime->pointer()[0][0]));
        diis->WriteErrorVector(&(ShalfGradShalf->pointer()[0][0]));
        diis->Extrapolate(&(Fprime->pointer()[0][0]));

        // Diagonalize F' to obtain C'
        Fprime->diagonalize(Fevec,Feval,ascending);

        // Find C = S^(-1/2)C'
        Ca_->gemm(false,false,1.0,Shalf,Fevec,0.0);

        outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",iter,energy_,dele,gnorm); 

        iter++;
        if ( iter > maxiter ) break;

    }while(fabs(dele) > e_convergence || gnorm > d_convergence );

    if ( iter > maxiter ) {
        throw PsiException("Maximum number of iterations exceeded!",__FILE__,__LINE__);
    }

    outfile->Printf("\n");
    outfile->Printf("    SCF iterations converged!\n");
    outfile->Printf("\n");

    outfile->Printf("    * Polaritonic RHF total energy: %20.12lf\n",energy_);

    // print cavity properties
    print_cavity_properties_ = true;
    build_cavity_hamiltonian();
    print_cavity_properties_ = false;
    
    Process::environment.globals["SCF TOTAL ENERGY"] = energy_;

    // print orbital energies
    Feval->print();

    return energy_;

}

} // End namespaces
