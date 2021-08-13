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

// for dft
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libscf_solver/hf.h"

// diis solver
#include <misc/diis.h>

#include "uks.h"

namespace hilbert{ 

PolaritonicUKS::PolaritonicUKS(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_):
    PolaritonicHF(reference_wavefunction,options_) {
    common_init();
}

PolaritonicUKS::~PolaritonicUKS() {
}

void PolaritonicUKS::common_init() {

    outfile->Printf("\n\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Polaritonic UKS                                  *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");

    // ensure scf_type df
    if ( options_.get_str("SCF_TYPE") != "DF" && options_.get_str("SCF_TYPE") != "CD" ) {
        throw PsiException("polaritonic uks only works with scf_type df for now",__FILE__,__LINE__);
    }

    // ensure running in c1 symmetry
    if ( reference_wavefunction_->nirrep() > 1 ) {
        throw PsiException("polaritonic uks only works with c1 symmetry for now.",__FILE__,__LINE__);
    }

    // SO-basis xc potential matrices
    Va_ = std::shared_ptr<Matrix>(new Matrix(nso_,nso_));
    Vb_ = std::shared_ptr<Matrix>(new Matrix(nso_,nso_));

    same_a_b_orbs_ = false;
    same_a_b_dens_ = false;

}

double PolaritonicUKS::compute_energy() {

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

    // determine the DFT functional and initialize the potential object
    psi::scf::HF* scfwfn = (psi::scf::HF*)reference_wavefunction_.get();
    std::shared_ptr<SuperFunctional> functional = scfwfn->functional();
    std::shared_ptr<VBase> potential = VBase::build_V(primary,functional,options_,(options_.get_str("REFERENCE") == "RKS" ? "RV" : "UV"));

    potential->initialize();

    // print the ks information
    potential->print_header();

    // JK object
    std::shared_ptr<JK> jk;

    int nQ = 0;
    bool is_x_lrc = false;
    if ( options_.get_str("SCF_TYPE") == "DF" ) {

        // get auxiliary basis:
        std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

        // total number of auxiliary basis functions
        nQ = auxiliary->nbf();

        std::shared_ptr<DiskDFJK> myjk = (std::shared_ptr<DiskDFJK>)(new DiskDFJK(primary,auxiliary));

        // memory for jk (say, 80% of what is available)
        myjk->set_memory(0.8 * memory_);

        // integral cutoff
        myjk->set_cutoff(options_.get_double("INTS_TOLERANCE"));

        // Do J/K/wK?
        is_x_lrc  = functional->is_x_lrc();
        //if ( options_["IP_FITTING"].has_changed() ) {
        //    if ( options_.get_bool("IP_FITTING") ) {
        //        is_x_lrc = true;
        //    }
        //}
        double x_omega = functional->x_omega();
        if ( options_["DFT_OMEGA"].has_changed() ) {
            x_omega = options_.get_double("DFT_OMEGA");
        }

        myjk->set_do_J(true);
        myjk->set_do_K(functional->is_x_hybrid());
        myjk->set_do_wK(is_x_lrc);
        myjk->set_omega(x_omega);

        myjk->initialize();

        jk = myjk;

    }else if ( options_.get_str("SCF_TYPE") == "CD" ) {

        std::shared_ptr<CDJK> myjk = (std::shared_ptr<CDJK>)(new CDJK(primary,options_.get_double("CHOLESKY_TOLERANCE")));

        // memory for jk (say, 80% of what is available)
        myjk->set_memory(0.8 * memory_);

        // integral cutoff
        myjk->set_cutoff(options_.get_double("INTS_TOLERANCE"));

        // Do J/K/wK?
        is_x_lrc  = functional->is_x_lrc();
        //if ( options_["IP_FITTING"].has_changed() ) {
        //    if ( options_.get_bool("IP_FITTING") ) {
        //        is_x_lrc = true;
        //    }
        //}
        double x_omega = functional->x_omega();
        if ( options_["DFT_OMEGA"].has_changed() ) {
            x_omega = options_.get_double("DFT_OMEGA");
        }

        myjk->set_do_J(true);
        myjk->set_do_K(functional->is_x_hybrid());
        myjk->set_do_wK(is_x_lrc);
        myjk->set_omega(x_omega);

        myjk->initialize();

        jk = myjk;

    }

    // grab some input options_
    double e_convergence = options_.get_double("E_CONVERGENCE");
    double d_convergence = options_.get_double("D_CONVERGENCE");
    int maxiter          = options_.get_int("MAXITER");

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso_);
    if ( options_.get_str("SCF_TYPE") == "DF" ) {
        outfile->Printf("    No. auxiliary basis functions:  %5i\n",nQ);
    }else if ( options_.get_str("SCF_TYPE") == "CD" ) {
        outfile->Printf("    cholesky_tolerance:             %5le\n",options_.get_double("CHOLESKY_TOLERANCE"));
    }
    outfile->Printf("    No. alpha electrons:            %5i\n",nalpha_);
    outfile->Printf("    No. beta electrons:             %5i\n",nbeta_);
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
    std::shared_ptr<Matrix> Fevec_a ( new Matrix(nso_,nso_) );
    std::shared_ptr<Matrix> Fevec_b ( new Matrix(nso_,nso_) );
    std::shared_ptr<Matrix> Fprime_a ( new Matrix(Fa_) );
    std::shared_ptr<Matrix> Fprime_b ( new Matrix(Fb_) );

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

    energy_  = enuc_ + average_electric_dipole_self_energy_;
    energy_ += 0.5 * Da_->vector_dot(h);
    energy_ += 0.5 * Db_->vector_dot(h);
    energy_ += 0.5 * Da_->vector_dot(Fa_);
    energy_ += 0.5 * Db_->vector_dot(Fb_);

    // SCF iterations

    double e_last    = 0.0;
    double dele      = 0.0;
    double gnorm_a   = 0.0;
    double gnorm_b   = 0.0;

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

    std::shared_ptr<DIIS> diis (new DIIS(2*nso_*nso_));

    int iter = 0;
    do {

        e_last = energy_;

        // grab occupied orbitals (the first nalpha)
        std::shared_ptr<Matrix> myCa (new Matrix(Ca_) );
        myCa->zero();

        // grab occupied orbitals (the first nbeta)
        std::shared_ptr<Matrix> myCb (new Matrix(Cb_) );
        myCb->zero();

        for (int mu = 0; mu < nso_; mu++) {
            for (int i = 0; i < nalpha_; i++) {
                myCa->pointer()[mu][i] = Ca_->pointer()[mu][i];
            }
            for (int i = 0; i < nbeta_; i++) {
                myCb->pointer()[mu][i] = Cb_->pointer()[mu][i];
            }
        }

        // push occupied orbitals onto JK object
        std::vector< std::shared_ptr<Matrix> >& C_left  = jk->C_left();
        C_left.clear();
        C_left.push_back(myCa);
        C_left.push_back(myCb);

        // form J/K
        jk->compute();

        // form Fa = h + Ja + Jb - Ka
        Fa_->copy(jk->J()[0]);
        Fa_->add(jk->J()[1]);

        Fb_->copy(jk->J()[0]);
        Fb_->add(jk->J()[1]);

        // Construct density from C
        C_DGEMM('n','t',nso_,nso_,nalpha_,1.0,&(Ca_->pointer()[0][0]),nso_,&(Ca_->pointer()[0][0]),nso_,0.0,&(Da_->pointer()[0][0]),nso_);
        C_DGEMM('n','t',nso_,nso_,nbeta_, 1.0,&(Cb_->pointer()[0][0]),nso_,&(Cb_->pointer()[0][0]),nso_,0.0,&(Db_->pointer()[0][0]),nso_);

        if (functional->needs_xc()) {

            // set a/b densities in potential object
            potential->set_D({Da_, Db_});

            // evaluate a/b potentials
            potential->compute_V({Va_,Vb_});

            // form Fa/b = h + Ja + Jb + Va/b
            Fa_->add(Va_);
            Fb_->add(Vb_);

        }

        // exact exchange?
        if (functional->is_x_hybrid()) {
            // form F = h + 2*J + V - alpha K
            double alpha = functional->x_alpha();
            Fa_->axpy(-alpha,jk->K()[0]);
            Fb_->axpy(-alpha,jk->K()[1]);
        }

        // LRC functional?
        if (is_x_lrc) {
            // form Fa/b = h + Ja + Jb + Va/b - alpha Ka/b - beta wKa/b
            double beta = 1.0 - functional->x_alpha();
            Fa_->axpy(-beta,jk->wK()[0]);
            Fb_->axpy(-beta,jk->wK()[1]);
        }

        std::shared_ptr<Matrix> oei (new Matrix(h));
        std::shared_ptr<Matrix> dipole_Ja (new Matrix(nso_,nso_));
        std::shared_ptr<Matrix> dipole_Jb (new Matrix(nso_,nso_));
        std::shared_ptr<Matrix> dipole_Ka (new Matrix(nso_,nso_));
        std::shared_ptr<Matrix> dipole_Kb (new Matrix(nso_,nso_));

        dipole_Ja->zero();
        dipole_Jb->zero();
        dipole_Ka->zero();
        dipole_Kb->zero();

        if ( n_photon_states_ > 1 ) {

            update_cavity_terms();

            // dipole self energy:

            // e-n term 
            oei->axpy(1.0,scaled_e_n_dipole_squared_);

            // one-electron part of e-e term 
            oei->axpy(-1.0,quadrupole_scaled_sum_);

            // two-electron part of e-e term (J)
            double scaled_mu_a = Da_->vector_dot(dipole_scaled_sum_);
            double scaled_mu_b = Db_->vector_dot(dipole_scaled_sum_);

            dipole_Ja->axpy(scaled_mu_a,dipole_scaled_sum_);
            dipole_Jb->axpy(scaled_mu_b,dipole_scaled_sum_);

            // two-electron part of e-e term (K)

            // Kpq += mu_pr * mu_qs * Drs
            double ** dp  = dipole_scaled_sum_->pointer();
            double ** dap = Da_->pointer();
            double ** dbp = Db_->pointer();
            double ** kap = dipole_Ka->pointer();
            double ** kbp = dipole_Kb->pointer();

            std::shared_ptr<Matrix> tmp (new Matrix(nso_,nso_));
            double ** tp = tmp->pointer();

            C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dp[0][0]),nso_,&(dap[0][0]),nso_,0.0,&(tp[0][0]),nso_);
            C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dp[0][0]),nso_,0.0,&(kap[0][0]),nso_);

            C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dp[0][0]),nso_,&(dbp[0][0]),nso_,0.0,&(tp[0][0]),nso_);
            C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dp[0][0]),nso_,0.0,&(kbp[0][0]),nso_);

            Fa_->add(dipole_Ja);
            Fa_->add(dipole_Jb);
            Fa_->subtract(dipole_Ka);

            Fb_->add(dipole_Ja);
            Fb_->add(dipole_Jb);
            Fb_->subtract(dipole_Kb);
        }

        Fa_->add(oei);
        Fb_->add(oei);

        // evaluate the current energy, E = D(H+F) + Enuc
        energy_  = enuc_ + average_electric_dipole_self_energy_;

        energy_ += Da_->vector_dot(oei);
        energy_ += Db_->vector_dot(oei);

        energy_ += 0.5 * Da_->vector_dot(jk->J()[0]);
        energy_ += 0.5 * Da_->vector_dot(jk->J()[1]);

        energy_ += 0.5 * Db_->vector_dot(jk->J()[0]);
        energy_ += 0.5 * Db_->vector_dot(jk->J()[1]);

        if (functional->is_x_hybrid()) {
            double alpha = functional->x_alpha();
            energy_ -= 0.5 * alpha * Da_->vector_dot(jk->K()[0]);
            energy_ -= 0.5 * alpha * Db_->vector_dot(jk->K()[1]);
        }

        energy_ += 0.5 * Da_->vector_dot(dipole_Ja);
        energy_ += 0.5 * Da_->vector_dot(dipole_Jb);
        energy_ -= 0.5 * Da_->vector_dot(dipole_Ka);

        energy_ += 0.5 * Db_->vector_dot(dipole_Ja);
        energy_ += 0.5 * Db_->vector_dot(dipole_Jb);
        energy_ -= 0.5 * Db_->vector_dot(dipole_Kb);

        if (is_x_lrc) {
            double beta = 1.0 - functional->x_alpha();
            energy_ -= 0.5 * beta * Da_->vector_dot(jk->wK()[0]);
            energy_ -= 0.5 * beta * Db_->vector_dot(jk->wK()[1]);
        }

        double exchange_correlation_energy = 0.0;
        if (functional->needs_xc()) {
            energy_ += potential->quadrature_values()["FUNCTIONAL"];
        }

        // dele
        dele = energy_ - e_last;

        // form F' = ST^(-1/2) F S^(-1/2)
        Fprime_a->copy(Fa_);
        Fprime_a->transform(Shalf);

        Fprime_b->copy(Fb_);
        Fprime_b->transform(Shalf);

        // The error vector in DIIS for SCF is defined as 
        // the orbital gradient, in the orthonormal basis:
        // 
        // ST^{-1/2} [FDS - SDF] S^{-1/2}

        std::shared_ptr<Matrix> grad_a = OrbitalGradient(Da_,Fa_,Shalf);
        std::shared_ptr<Matrix> grad_b = OrbitalGradient(Db_,Fb_,Shalf);

        // We will use the RMS of the orbital gradient 
        // to monitor convergence.
        gnorm_a = grad_a->rms();
        gnorm_b = grad_b->rms();

        // DIIS extrapolation
        diis->WriteVector(&(Fprime_a->pointer()[0][0]),&(Fprime_b->pointer()[0][0]));
        diis->WriteErrorVector(&(grad_a->pointer()[0][0]),&(grad_b->pointer()[0][0]));
        diis->Extrapolate(&(Fprime_a->pointer()[0][0]),&(Fprime_b->pointer()[0][0]));

        // Diagonalize F' to obtain C'
        Fprime_a->diagonalize(Fevec_a,epsilon_a_,ascending);
        Fprime_b->diagonalize(Fevec_b,epsilon_b_,ascending);

        // Find C = S^(-1/2)C'
        Ca_->gemm(false,false,1.0,Shalf,Fevec_a,0.0);
        Cb_->gemm(false,false,1.0,Shalf,Fevec_b,0.0);

        outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",iter,energy_,dele,0.5 * (gnorm_a + gnorm_b)); 

        iter++;
        if ( iter > maxiter ) break;

    }while(fabs(dele) > e_convergence || 0.5 * (gnorm_a + gnorm_b) > d_convergence );

    if ( iter > maxiter ) {
        throw PsiException("Maximum number of iterations exceeded!",__FILE__,__LINE__);
    }

    outfile->Printf("\n");
    outfile->Printf("    SCF iterations converged!\n");
    outfile->Printf("\n");

    outfile->Printf("    * Polaritonic UKS total energy: %20.12lf\n",energy_);

    Process::environment.globals["SCF TOTAL ENERGY"] = energy_;

    Process::environment.globals["CURRENT ENERGY"] = energy_;

    // update cavity terms once more
    if ( n_photon_states_ > 1 ) {
        update_cavity_terms();
    } 

    // print orbital energies
    epsilon_a_->print();
    epsilon_b_->print();

    return energy_;

}

} // End namespaces
