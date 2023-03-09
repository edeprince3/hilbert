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

#include "rks.h"

namespace hilbert{ 

PolaritonicRKS::PolaritonicRKS(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_):
    PolaritonicHF(reference_wavefunction,options_) {
    common_init();
}

PolaritonicRKS::~PolaritonicRKS() {
}

void PolaritonicRKS::common_init() {

    outfile->Printf("\n\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Polaritonic RKS                                  *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");

    // ensure scf_type df
    if ( options_.get_str("SCF_TYPE") != "DF" && options_.get_str("SCF_TYPE") != "CD" ) {
        throw PsiException("polaritonic rks only works with scf_type df for now",__FILE__,__LINE__);
    }

    // ensure running in c1 symmetry
    if ( reference_wavefunction_->nirrep() > 1 ) {
        throw PsiException("polaritonic rks only works with c1 symmetry for now.",__FILE__,__LINE__);
    }

    // SO-basis xc potential matrices
    Va_ = std::shared_ptr<Matrix>(new Matrix(nso_,nso_));

    same_a_b_orbs_ = true;
    same_a_b_dens_ = true;

    n_photon_states_ = options_.get_int("N_PHOTON_STATES");

    // cavity Hamiltonian in basis of photon number states
    HCavity_x_               = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    HCavity_y_               = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    HCavity_z_               = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));

    // molecule->cavity interactoin Hamiltonian in basis of photon number states
    HCavityInteraction_x_    = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    HCavityInteraction_y_    = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    HCavityInteraction_z_    = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));

    // total cavity Hamiltonian in basis of photon number states
    HCavityTotal_x_          = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    HCavityTotal_y_          = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    HCavityTotal_z_          = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));

    // cavity dipole operator in basis of photon number states
    CavityDipole_x_          = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    CavityDipole_y_          = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));
    CavityDipole_z_          = (std::shared_ptr<Matrix>)(new Matrix(n_photon_states_,n_photon_states_));

    use_coherent_state_basis_ = options_.get_bool("USE_COHERENT_STATE_BASIS");

}

double PolaritonicRKS::compute_energy() {

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
    //std::shared_ptr<JK> jk;

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

        jk_ = myjk;

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

        jk_ = myjk;

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
    std::shared_ptr<Matrix> Fprime_a ( new Matrix(Fa_) );

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

    if ( !options_.get_bool("QED_USE_RELAXED_ORBITALS") && use_coherent_state_basis_ ) {
        throw PsiException("QED_USE_RELAXED_ORBITALS false is not compatible with USE_COHERENT_STATE_BASIS true",__FILE__,__LINE__);
    }

    energy_  = enuc_;
    if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {
        energy_  += average_electric_dipole_self_energy_;
    }
    energy_ += 0.5 * Da_->vector_dot(h);
    energy_ += 0.5 * Da_->vector_dot(h);
    energy_ += 0.5 * Da_->vector_dot(Fa_);
    energy_ += 0.5 * Da_->vector_dot(Fa_);

    // SCF iterations

    double e_last    = 0.0;
    double dele      = 0.0;
    double gnorm_a   = 0.0;

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

        // w <b*b>
        double cavity_energy = 0.0;

        e_last = energy_;

        // grab occupied orbitals (the first nalpha)
        std::shared_ptr<Matrix> myCa (new Matrix(Ca_) );
        myCa->zero();

        for (int mu = 0; mu < nso_; mu++) {
            for (int i = 0; i < nalpha_; i++) {
                myCa->pointer()[mu][i] = Ca_->pointer()[mu][i];
            }
        }

        // push occupied orbitals onto JK object
        std::vector< std::shared_ptr<Matrix> >& C_left  = jk_->C_left();
        C_left.clear();
        C_left.push_back(myCa);

        // form J/K
        jk_->compute();

        // form Fa = h + Ja + Jb - Ka
        Fa_->copy(jk_->J()[0]);
        Fa_->add(jk_->J()[0]);

        // Construct density from C
        C_DGEMM('n','t',nso_,nso_,nalpha_,1.0,&(Ca_->pointer()[0][0]),nso_,&(Ca_->pointer()[0][0]),nso_,0.0,&(Da_->pointer()[0][0]),nso_);

        if (functional->needs_xc()) {

            // set a/b densities in potential object
            potential->set_D({Da_});

            // evaluate a/b potentials
            potential->compute_V({Va_});

            // form Fa/b = h + Ja + Jb + Va/b
            Fa_->add(Va_);

        }

        // exact exchange?
        if (functional->is_x_hybrid()) {
            // form F = h + 2*J + V - alpha K
            double alpha = functional->x_alpha();
            Fa_->axpy(-alpha,jk_->K()[0]);
        }

        // LRC functional?
        if (is_x_lrc) {
            // form Fa/b = h + Ja + Jb + Va/b - alpha Ka/b - beta wKa/b
            double beta = 1.0 - functional->x_alpha();
            Fa_->axpy(-beta,jk_->wK()[0]);
        }

        std::shared_ptr<Matrix> oei (new Matrix(h));
        std::shared_ptr<Matrix> dipole_Ja (new Matrix(nso_,nso_));
        std::shared_ptr<Matrix> dipole_Jb (new Matrix(nso_,nso_));
        std::shared_ptr<Matrix> dipole_Ka (new Matrix(nso_,nso_));
        std::shared_ptr<Matrix> dipole_Kb (new Matrix(nso_,nso_));

        dipole_Ja->zero();
        dipole_Ka->zero();

        if ( n_photon_states_ > 1 ) {

            update_cavity_terms();

            if ( !use_coherent_state_basis_ ) {

                // build and diagonalize cavity hamiltonian (for canonical basis only)
                cavity_energy = build_cavity_hamiltonian();

                // photon -> electron contribution (for canonical basis only)
                // 
                // -(w/2)^1/2 lambda.mu_e < mu_p >
                // 
                std::shared_ptr<Matrix> V = (std::shared_ptr<Matrix>)(new Matrix(dipole_scaled_sum_));
                V->scale(-sqrt(0.5 * cavity_frequency_[2])*CavityDipole_z_->pointer()[0][0]);
                oei->add(V);
            }

            // dipole self energy ... ignore if following QED-TDDFT outlined in J. Chem. Phys. 155, 064107 (2021)

            if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {

                // e-n term 
                oei->axpy(1.0,scaled_e_n_dipole_squared_);

                // one-electron part of e-e term 
                oei->axpy(-1.0,quadrupole_scaled_sum_);

                // two-electron part of e-e term (J)
                double scaled_mu_a = Da_->vector_dot(dipole_scaled_sum_);

                dipole_Ja->axpy(scaled_mu_a,dipole_scaled_sum_);

                // two-electron part of e-e term (K)

                // Kpq += mu_pr * mu_qs * Drs
                double ** dp  = dipole_scaled_sum_->pointer();
                double ** dap = Da_->pointer();
                double ** kap = dipole_Ka->pointer();

                std::shared_ptr<Matrix> tmp (new Matrix(nso_,nso_));
                double ** tp = tmp->pointer();

                C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dp[0][0]),nso_,&(dap[0][0]),nso_,0.0,&(tp[0][0]),nso_);
                C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dp[0][0]),nso_,0.0,&(kap[0][0]),nso_);

                Fa_->add(dipole_Ja);
                Fa_->add(dipole_Ja);
                Fa_->subtract(dipole_Ka);

            }

        }

        Fa_->add(oei);

        // evaluate the current energy, E = D(H+F) + Enuc
        // 
        // note that, in canonical basis, average_electric_dipole_self_energy_ 
        // only contains nuclear contributions
        // 
        energy_  = enuc_;
        if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {
            energy_  += average_electric_dipole_self_energy_;
        }

        energy_ += Da_->vector_dot(oei);
        energy_ += Da_->vector_dot(oei);

        energy_ += 0.5 * Da_->vector_dot(jk_->J()[0]);
        energy_ += 0.5 * Da_->vector_dot(jk_->J()[0]);

        energy_ += 0.5 * Da_->vector_dot(jk_->J()[0]);
        energy_ += 0.5 * Da_->vector_dot(jk_->J()[0]);

        if (functional->is_x_hybrid()) {
            double alpha = functional->x_alpha();
            energy_ -= 0.5 * alpha * Da_->vector_dot(jk_->K()[0]);
            energy_ -= 0.5 * alpha * Da_->vector_dot(jk_->K()[0]);
        }

        // dipole self energy contributions ... ignore if following QED-TDDFT outlined in J. Chem. Phys. 155, 064107 (2021)
        if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {
            energy_ += 0.5 * Da_->vector_dot(dipole_Ja);
            energy_ += 0.5 * Da_->vector_dot(dipole_Ja);
            energy_ -= 0.5 * Da_->vector_dot(dipole_Ka);
    
            energy_ += 0.5 * Da_->vector_dot(dipole_Ja);
            energy_ += 0.5 * Da_->vector_dot(dipole_Ja);
            energy_ -= 0.5 * Da_->vector_dot(dipole_Ka);
        }
    
        if (is_x_lrc) {
            double beta = 1.0 - functional->x_alpha();
            energy_ -= 0.5 * beta * Da_->vector_dot(jk_->wK()[0]);
            energy_ -= 0.5 * beta * Da_->vector_dot(jk_->wK()[0]);
        }

        double exchange_correlation_energy = 0.0;
        if (functional->needs_xc()) {
            energy_ += potential->quadrature_values()["FUNCTIONAL"];
        }

        // if using canonical basis, add additional terms. ignore if following QED-TDDFT outlined in J. Chem. Phys. 155, 064107 (2021)
        if ( !use_coherent_state_basis_ && options_.get_bool("QED_USE_RELAXED_ORBITALS")) {

            // w <b*b> 
            energy_ += cavity_energy; 

            // -(w/2)^1/2 lambda.mu_n <b* + b> 
            double lambda_z = cavity_coupling_strength_[2] * sqrt(2.0 * cavity_frequency_[2]);
            energy_ += -sqrt(0.5 * cavity_frequency_[2]) * lambda_z * nuc_dip_z_ * CavityDipole_z_->pointer()[0][0];

        }

        // dele
        dele = energy_ - e_last;

        // form F' = ST^(-1/2) F S^(-1/2)
        Fprime_a->copy(Fa_);
        Fprime_a->transform(Shalf);

        // The error vector in DIIS for SCF is defined as 
        // the orbital gradient, in the orthonormal basis:
        // 
        // ST^{-1/2} [FDS - SDF] S^{-1/2}

        std::shared_ptr<Matrix> grad_a = OrbitalGradient(Da_,Fa_,Shalf);

        // We will use the RMS of the orbital gradient 
        // to monitor convergence.
        gnorm_a = grad_a->rms();

        // DIIS extrapolation
        diis->WriteVector(&(Fprime_a->pointer()[0][0]));
        diis->WriteErrorVector(&(grad_a->pointer()[0][0]));
        diis->Extrapolate(&(Fprime_a->pointer()[0][0]));

        // Diagonalize F' to obtain C'
        Fprime_a->diagonalize(Fevec_a,epsilon_a_,ascending);

        // Find C = S^(-1/2)C'
        Ca_->gemm(false,false,1.0,Shalf,Fevec_a,0.0);

        outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",iter,energy_,dele,gnorm_a); 

        iter++;
        if ( iter > maxiter ) break;

    }while(fabs(dele) > e_convergence || gnorm_a > d_convergence );

    if ( iter > maxiter ) {
        throw PsiException("Maximum number of iterations exceeded!",__FILE__,__LINE__);
    }

    outfile->Printf("\n");
    outfile->Printf("    SCF iterations converged!\n");
    outfile->Printf("\n");

    // evaluate dipole self energy
    if ( n_photon_states_ > 1 ) {
        evaluate_dipole_self_energy();
    }

    evaluate_dipole_variance();

    outfile->Printf("    * Polaritonic RKS total energy: %20.12lf\n",energy_);
    outfile->Printf("\n");

    Process::environment.globals["SCF TOTAL ENERGY"] = energy_;

    Process::environment.globals["CURRENT ENERGY"] = energy_;

    // update cavity terms once more
    if ( n_photon_states_ > 1 ) {
        update_cavity_terms();
    } 

    // print orbital energies
    epsilon_a_->print();

    return energy_;

}

/* 
 * build_cavity_hamiltonian():
 *
 * Build the cavity Hamiltonian, the molecule->cavity interaction
 * Hamiltonian, and the cavity dipole moment operator in the
 * basis of photon number states.  Transform all three operators
 * to the basis that diagonalizes the total cavity Hamiltonian
 * (the sum of the cavity Hamiltonian and the molecule->cavity
 * Hamiltonian) and returns its lowest eigenvalue (in z-direction)
*/
double PolaritonicRKS::build_cavity_hamiltonian(){

    // First, build the cavity dipole moment operator 
    // in the basis of photon number states ( n_photon_states_ )

    int nS = n_photon_states_;

    CavityDipole_x_->zero();
    CavityDipole_y_->zero();
    CavityDipole_z_->zero();

    for (int A=0; A<nS-1; A++){
        //CavityDipole_x_->pointer()[A+1][A] += cavity_frequency_[0] * cavity_coupling_strength_[0] * sqrt(A+1);
        //CavityDipole_x_->pointer()[A][A+1] += cavity_frequency_[0] * cavity_coupling_strength_[0] * sqrt(A+1);
        //CavityDipole_y_->pointer()[A+1][A] += cavity_frequency_[1] * cavity_coupling_strength_[1] * sqrt(A+1);
        //CavityDipole_y_->pointer()[A][A+1] += cavity_frequency_[1] * cavity_coupling_strength_[1] * sqrt(A+1);
        //CavityDipole_z_->pointer()[A+1][A] += cavity_frequency_[2] * cavity_coupling_strength_[2] * sqrt(A+1);
        //CavityDipole_z_->pointer()[A][A+1] += cavity_frequency_[2] * cavity_coupling_strength_[2] * sqrt(A+1);
        //CavityDipole_x_->pointer()[A+1][A] += sqrt(A+1);
        //CavityDipole_x_->pointer()[A][A+1] += sqrt(A+1);
        //CavityDipole_y_->pointer()[A+1][A] += sqrt(A+1);
        //CavityDipole_y_->pointer()[A][A+1] += sqrt(A+1);
        CavityDipole_z_->pointer()[A+1][A] += sqrt(A+1);
        CavityDipole_z_->pointer()[A][A+1] += sqrt(A+1);
    }

    // Build the cavity Hamiltonian operator in the basis of 
    // photon number states ( n_photon_states_ )

    HCavity_x_->zero();
    HCavity_y_->zero();
    HCavity_z_->zero();

    for (int A=0; A<nS; A++){
        HCavity_x_->pointer()[A][A] = A*cavity_frequency_[0];
        HCavity_y_->pointer()[A][A] = A*cavity_frequency_[1];
        HCavity_z_->pointer()[A][A] = A*cavity_frequency_[2];
    }

    // Evaluate the electronic contribute to the molecule's dipole moment

    double e_dip_x = C_DDOT(nso_*nso_,&(Da_->pointer())[0][0],1,&(dipole_[0]->pointer())[0][0],1); 
    double e_dip_y = C_DDOT(nso_*nso_,&(Da_->pointer())[0][0],1,&(dipole_[1]->pointer())[0][0],1); 
    double e_dip_z = C_DDOT(nso_*nso_,&(Da_->pointer())[0][0],1,&(dipole_[2]->pointer())[0][0],1);
    if ( same_a_b_dens_ ) {
        e_dip_x *= 2.0;
        e_dip_y *= 2.0;
        e_dip_z *= 2.0;
    }else {
        e_dip_x += C_DDOT(nso_*nso_,&(Db_->pointer())[0][0],1,&(dipole_[0]->pointer())[0][0],1); 
        e_dip_y += C_DDOT(nso_*nso_,&(Db_->pointer())[0][0],1,&(dipole_[1]->pointer())[0][0],1); 
        e_dip_z += C_DDOT(nso_*nso_,&(Db_->pointer())[0][0],1,&(dipole_[2]->pointer())[0][0],1);
    }

    // Build the molecule->cavity interaction Hamiltonian operator 
    // in the basis of photon number states ( n_photon_states_ )

    HCavityInteraction_x_->zero();
    HCavityInteraction_y_->zero();
    HCavityInteraction_z_->zero();

    for (int A = 0; A<nS; A++) {
        if (A < nS - 1) {

            HCavityInteraction_x_->pointer()[A][A+1] = -(e_dip_x + nuc_dip_x_) * cavity_frequency_[0] * cavity_coupling_strength_[0] * sqrt(A+1);
            HCavityInteraction_y_->pointer()[A][A+1] = -(e_dip_y + nuc_dip_y_) * cavity_frequency_[1] * cavity_coupling_strength_[1] * sqrt(A+1);
            HCavityInteraction_z_->pointer()[A][A+1] = -(e_dip_z + nuc_dip_z_) * cavity_frequency_[2] * cavity_coupling_strength_[2] * sqrt(A+1);

        }

        if (A > 0) {

            HCavityInteraction_x_->pointer()[A][A-1] = -(e_dip_x + nuc_dip_x_) * cavity_frequency_[0] * cavity_coupling_strength_[0] * sqrt(A);
            HCavityInteraction_y_->pointer()[A][A-1] = -(e_dip_y + nuc_dip_y_) * cavity_frequency_[1] * cavity_coupling_strength_[1] * sqrt(A);
            HCavityInteraction_z_->pointer()[A][A-1] = -(e_dip_z + nuc_dip_z_) * cavity_frequency_[2] * cavity_coupling_strength_[2] * sqrt(A);
        }
    }

    // Build the total cavity Hamiltonian

    HCavityTotal_x_->copy(HCavity_x_);
    HCavityTotal_y_->copy(HCavity_y_);
    HCavityTotal_z_->copy(HCavity_z_);
    HCavityTotal_x_->add(HCavityInteraction_x_);
    HCavityTotal_y_->add(HCavityInteraction_y_);
    HCavityTotal_z_->add(HCavityInteraction_z_);

    // Diagonalize total cavity Hamiltonian and transform cavity
    // Hamiltonian, molecule->cavity interaction Hamiltonian, and
    // cavity dipole moment operator to the basis that diagonalizes
    // the total cavity Hamiltonian
    std::shared_ptr<Matrix> eigvec (new Matrix(n_photon_states_,n_photon_states_));
    std::shared_ptr<Vector> eigval (new Vector(n_photon_states_));

    // x
    HCavityTotal_x_->diagonalize(eigvec, eigval);
    HCavity_x_->transform(eigvec);
    HCavityInteraction_x_->transform(eigvec);
    CavityDipole_x_->transform(eigvec);
    double ex = eigval->pointer()[0];

    // y
    HCavityTotal_y_->diagonalize(eigvec, eigval);
    HCavity_y_->transform(eigvec);
    HCavityInteraction_y_->transform(eigvec);
    CavityDipole_y_->transform(eigvec);
    double ey = eigval->pointer()[0];

    // z
    HCavityTotal_z_->diagonalize(eigvec, eigval);
    HCavity_z_->transform(eigvec);
    HCavityInteraction_z_->transform(eigvec);
    CavityDipole_z_->transform(eigvec);
    double ez = eigval->pointer()[0];

    if ( print_cavity_properties_ ) {
        outfile->Printf("\n");
        outfile->Printf("    ==> Cavity Properties <==\n");
        outfile->Printf("\n");

        outfile->Printf("        Energy x: %20.12lf    Dipole x: %20.12lf\n",ex,CavityDipole_x_->pointer()[0][0]);
        outfile->Printf("        Energy y: %20.12lf    Dipole y: %20.12lf\n",ey,CavityDipole_y_->pointer()[0][0]);
        outfile->Printf("        Energy z: %20.12lf    Dipole z: %20.12lf\n",ez,CavityDipole_z_->pointer()[0][0]);
        outfile->Printf("\n");
    }

    return ez;

}

} // End namespaces
