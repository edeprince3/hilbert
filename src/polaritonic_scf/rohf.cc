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

#include "rohf.h"

namespace hilbert{ 

PolaritonicROHF::PolaritonicROHF(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_):
    PolaritonicHF(reference_wavefunction,options_) {
    common_init();
}

PolaritonicROHF::~PolaritonicROHF() {
}

void PolaritonicROHF::common_init() {

    outfile->Printf("\n\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Polaritonic ROHF                                 *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");

    same_a_b_orbs_ = true;
    same_a_b_dens_ = false;

    // ensure scf_type df
    if ( options_.get_str("SCF_TYPE") != "DF" && options_.get_str("SCF_TYPE") != "CD" ) {
        throw PsiException("polaritonic rohf only works with scf_type df for now",__FILE__,__LINE__);
    }

    // ensure running in c1 symmetry
    if ( reference_wavefunction_->nirrep() > 1 ) {
        throw PsiException("polaritonic rohf only works with c1 symmetry for now.",__FILE__,__LINE__);
    }

    // this code only works in the coheren-state basis
    if ( !options_.get_bool("USE_COHERENT_STATE_BASIS") ) {
        // but there is only a problem if we are relaxing the orbitals:
        if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {
            throw PsiException("polaritonic uhf can only relax orbitals within the coherent-state basis", __FILE__, __LINE__);
        }
    }

}

double PolaritonicROHF::compute_energy() {

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

    // JK object
    std::shared_ptr<JK> jk;

    int nQ = 0;
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

        // Do J/K, Not wK 
        myjk->set_do_J(true);
        myjk->set_do_K(true);
        myjk->set_do_wK(false);

        myjk->initialize();

        jk = myjk;

    }else if ( options_.get_str("SCF_TYPE") == "CD" ) {

        std::shared_ptr<CDJK> myjk = (std::shared_ptr<CDJK>)(new CDJK(primary,options_.get_double("CHOLESKY_TOLERANCE")));

        // memory for jk (say, 80% of what is available)
        myjk->set_memory(0.8 * memory_);

        // integral cutoff
        myjk->set_cutoff(options_.get_double("INTS_TOLERANCE"));

        // Do J/K, Not wK 
        myjk->set_do_J(true);
        myjk->set_do_K(true);
        myjk->set_do_wK(false);

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

    std::shared_ptr<Matrix> Shalf ( new Matrix(nso_,nso_) );
    for (int mu = 0; mu < nso_; mu++) {
        Shalf->pointer()[mu][mu] = 1.0 / sqrt(Seval->pointer()[mu]);
    }

    // transform Seval back to nonorthogonal basis
    Shalf->back_transform(Sevec);

    // eigenfunctions of orthonormal SO basis effective fock matrix
    std::shared_ptr<Matrix> Fevec = (std::shared_ptr<Matrix>)(new Matrix(nso_,nso_)); 

    // allocate memory for F' and its eigenvectors and eigenvalues
    std::shared_ptr<Matrix> Fprime ( new Matrix(nso_,nso_) );

    // allocate memory for effective Fock matrix
    std::shared_ptr<Matrix> Feff ( new Matrix(nso_,nso_) );

    // guess (from existing reference)

    //if ( !options_.get_bool("QED_USE_RELAXED_ORBITALS") && !use_coherent_state_basis_ ) {
    //    throw PsiException("QED_USE_RELAXED_ORBITALS false is not compatible with USE_COHERENT_STATE_BASIS false",__FILE__,__LINE__);
    //}

    // form initial Feff
    form_Feff(Feff);

    // form F' = ST^(-1/2) F S^(-1/2), where F = whatever was passed in
    Fprime->copy(Feff);
    Fprime->transform(Shalf);
    Fprime->diagonalize(Fevec, epsilon_a_,ascending);

    // Find C = S^(-1/2)C'
    Ca_->gemm(false, false, 1.0, Shalf, Fevec, 0.0);

    // Construct density from C
    C_DGEMM('n','t',nso_,nso_,nalpha_,1.0,&(Ca_->pointer()[0][0]),nso_,&(Ca_->pointer()[0][0]),nso_,0.0,&(Da_->pointer()[0][0]),nso_);
    C_DGEMM('n','t',nso_,nso_,nbeta_,1.0,&(Ca_->pointer()[0][0]),nso_,&(Ca_->pointer()[0][0]),nso_,0.0,&(Db_->pointer()[0][0]),nso_);

    energy_  = enuc_;
    if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {
        energy_  += average_electric_dipole_self_energy_;
    }
    energy_ += 0.5 * Da_->vector_dot(h);
    energy_ += 0.5 * Db_->vector_dot(h);
    energy_ += 0.5 * Da_->vector_dot(Fa_);
    energy_ += 0.5 * Db_->vector_dot(Fb_);

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
                myCb->pointer()[mu][i] = Ca_->pointer()[mu][i];
            }
        }

        // push occupied orbitals onto JK object
        std::vector< std::shared_ptr<Matrix> >& C_left  = jk->C_left();
        C_left.clear();
        C_left.push_back(myCa);
        C_left.push_back(myCb);

        // form J/K
        jk->compute();

        std::shared_ptr<Matrix> oei (new Matrix(h));

        if ( n_photon_states_ > 1 ) {

            update_cavity_terms();

            // dipole self energy ... ignore if not relaxing the orbitals

            if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {

                // e-n term: 1/2 lambda^2 * 2 de ( dn - <d> )
                oei->axpy(1.0,scaled_e_n_dipole_squared_);

                // one-electron part of e-e term 
                oei->axpy(-1.0,quadrupole_scaled_sum_);

                // two-electron part of e-e term (J)
                double scaled_mu_a = Da_->vector_dot(dipole_scaled_sum_);
                double scaled_mu_b = Db_->vector_dot(dipole_scaled_sum_);

                jk->J()[0]->axpy(scaled_mu_a,dipole_scaled_sum_);
                jk->J()[1]->axpy(scaled_mu_b,dipole_scaled_sum_);

                //Fa_->axpy(scaled_mu_a,dipole_scaled_sum_);
                //Fa_->axpy(scaled_mu_b,dipole_scaled_sum_);

                //Fb_->axpy(scaled_mu_a,dipole_scaled_sum_);
                //Fb_->axpy(scaled_mu_b,dipole_scaled_sum_);

                // two-electron part of e-e term (K)

                // Kpq += mu_pr * mu_qs * Drs
                double ** dp  = dipole_scaled_sum_->pointer();
                double ** dap = Da_->pointer();
                double ** dbp = Db_->pointer();
                double ** kap = jk->K()[0]->pointer();
                double ** kbp = jk->K()[1]->pointer();

                std::shared_ptr<Matrix> tmp (new Matrix(nso_,nso_));
                double ** tp = tmp->pointer();

                C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dp[0][0]),nso_,&(dap[0][0]),nso_,0.0,&(tp[0][0]),nso_);
                C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dp[0][0]),nso_,1.0,&(kap[0][0]),nso_);

                C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(dp[0][0]),nso_,&(dbp[0][0]),nso_,0.0,&(tp[0][0]),nso_);
                C_DGEMM('n','t',nso_,nso_,nso_,1.0,&(tp[0][0]),nso_,&(dp[0][0]),nso_,1.0,&(kbp[0][0]),nso_);

            }
        }

        // form Fa = h + Ja + Jb - Ka
        Fa_->copy(oei);
        Fa_->add(jk->J()[0]);
        Fa_->add(jk->J()[1]);
        Fa_->subtract(jk->K()[0]);

        Fb_->copy(oei);
        Fb_->add(jk->J()[0]);
        Fb_->add(jk->J()[1]);
        Fb_->subtract(jk->K()[1]);

        // Construct density from C
        C_DGEMM('n','t',nso_,nso_,nalpha_,1.0,&(Ca_->pointer()[0][0]),nso_,&(Ca_->pointer()[0][0]),nso_,0.0,&(Da_->pointer()[0][0]),nso_);
        C_DGEMM('n','t',nso_,nso_,nbeta_,1.0,&(Ca_->pointer()[0][0]),nso_,&(Ca_->pointer()[0][0]),nso_,0.0,&(Db_->pointer()[0][0]),nso_);

        // evaluate the current energy, E = D(H+F) + Enuc
        energy_  = enuc_;
        if ( options_.get_bool("QED_USE_RELAXED_ORBITALS") ) {
            energy_  += average_electric_dipole_self_energy_;
        }

        energy_ += Da_->vector_dot(oei);
        energy_ += Db_->vector_dot(oei);

        energy_ += 0.5 * Da_->vector_dot(jk->J()[0]);
        energy_ += 0.5 * Da_->vector_dot(jk->J()[1]);
        energy_ -= 0.5 * Da_->vector_dot(jk->K()[0]);

        energy_ += 0.5 * Db_->vector_dot(jk->J()[0]);
        energy_ += 0.5 * Db_->vector_dot(jk->J()[1]);
        energy_ -= 0.5 * Db_->vector_dot(jk->K()[1]);

        // dele
        dele = energy_ - e_last;

        // form effective Fock matrix in non-orthogonal SO basis
        form_Feff(Feff);

        // form orbital gradient in orthogonal SO basis
        Feff->transform(Shalf);
        std::shared_ptr<Matrix> grad = OrbitalGradientROHF(Feff,Fevec);

        // We will use the RMS of the orbital gradient 
        // to monitor convergence.
        gnorm = grad->rms();

        // form F' = ST^(-1/2) F S^(-1/2) for DIIS
        Fprime->copy(Feff);
        //Fprime->transform(Shalf);

        // DIIS extrapolation
        diis->WriteVector(&(Fprime->pointer()[0][0]));
        diis->WriteErrorVector(&(grad->pointer()[0][0]));
        diis->Extrapolate(&(Fprime->pointer()[0][0]));

        // Diagonalize F' to obtain C'
        Fprime->diagonalize(Fevec,epsilon_a_,ascending);

        // Find C = S^(-1/2)C'
        Ca_->gemm(false,false,1.0,Shalf,Fevec,0.0);
        Cb_->gemm(false,false,1.0,Shalf,Fevec,0.0);

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

    // evaluate dipole self energy
    if ( n_photon_states_ > 1 ) {
        evaluate_dipole_self_energy();
    }

    evaluate_dipole_variance();

    outfile->Printf("    * Polaritonic ROHF total energy: %20.12lf\n",energy_);

    Process::environment.globals["SCF TOTAL ENERGY"] = energy_;

    Process::environment.globals["CURRENT ENERGY"] = energy_;

    // print orbital energies
    Fprime->copy(Fa_);
    Fprime->transform(Shalf);
    Fprime->diagonalize(Fevec,epsilon_a_,ascending);
    epsilon_a_->print();

    Fprime->copy(Fb_);
    Fprime->transform(Shalf);
    Fprime->diagonalize(Fevec,epsilon_b_,ascending);
    epsilon_b_->print();

    // update cavity terms once more
    if ( n_photon_states_ > 1 ) {
        update_cavity_terms();
    } 

    return energy_;

}

void PolaritonicROHF::form_Feff(std::shared_ptr<Matrix> Feff) {

    /*
     * Fo = open-shell fock matrix = 0.5 Fa
     * Fc = closed-shell fock matrix = 0.5 (Fa + Fb)
     *
     * Therefore
     *
     * 2(Fc-Fo) = Fb
     * 2Fo = Fa
     *
     * Form the effective Fock matrix, too
     * The effective Fock matrix has the following structure
     *          |  closed     open    virtual
     *  ----------------------------------------
     *  closed  |    Fc     2(Fc-Fo)    Fc
     *  open    | 2(Fc-Fo)     Fc      2Fo
     *  virtual |    Fc       2Fo       Fc
     */

    std::shared_ptr<Matrix> moFa (new Matrix(Fa_));
    std::shared_ptr<Matrix> moFb (new Matrix(Fb_)); 

    moFa->transform(Fa_,Ca_);
    moFb->transform(Fb_,Ca_); 

    Feff->copy(moFa);
    Feff->add(moFb);
    Feff->scale(0.5);

    for (int h = 0; h < nirrep_; ++h) {
        for (int i = doccpi_[h]; i < doccpi_[h] + soccpi_[h]; ++i) {
            // Set the open/closed portion
            for (int j = 0; j < doccpi_[h]; ++j) {
                double val = moFb->get(h, i, j);
                Feff->set(h, i, j, val);
                Feff->set(h, j, i, val);
            }
            // Set the open/virtual portion
            for (int j = doccpi_[h] + soccpi_[h]; j < nmopi_[h]; ++j) {
                double val = moFa->get(h, i, j);
                Feff->set(h, i, j, val);
                Feff->set(h, j, i, val);
            }
        }
    }

    // Form the orthogonalized SO basis moFeff matrix, for use in DIIS
    //std::shared_ptr<Matrix> diag_F_temp (new Matrix(Fa_));
    //diag_F_temp->gemm(false, false, 1.0, Ct_, Feff, 0.0);
    //Feff->gemm(false, true, 1.0, diag_F_temp, Ct_, 0.0);

    // Form the no-orthogonalized SO basis moFeff matrix, for use in DIIS
    //Feff->back_transform(Ct_);
    //Feff->back_transform(Shalf2_);
    Feff->back_transform(Ca_);
    Feff->back_transform(S_);

    //Feff->back_transform(Shalf2_);

    //Feff->back_transform(Ct_);
    //Feff->back_transform(Ct_);

    //Feff->transform(Ca_);
    //Feff->print();
    //exit(0);

}

std::shared_ptr<Matrix> PolaritonicROHF::OrbitalGradientROHF(std::shared_ptr<Matrix> F, std::shared_ptr<Matrix> C) {

    // Only the inact-act, inact-vir, and act-vir rotations are non-redundant
    Dimension dim_zero = Dimension(nirrep_, "Zero Dim");
    Dimension noccpi = doccpi_ + soccpi_;
    Dimension virpi = nmopi_ - doccpi_;
    Slice row_slice(dim_zero, noccpi);
    Slice col_slice(doccpi_, doccpi_ + virpi);

    std::shared_ptr<Matrix> tmp (new Matrix(nso_,nso_));
    tmp->copy(F);
    tmp->transform(C);

    std::shared_ptr<Matrix> MOgradient = tmp->get_block(row_slice, col_slice);

    // Zero out act-act part
    for (size_t h = 0; h < nirrep_; h++) {
        if (!soccpi_[h]) continue;

        for (size_t i = 0; i < soccpi_[h]; i++) {
            for (size_t j = 0; j < soccpi_[h]; j++) {
                MOgradient->set(h, i + doccpi_[h], j, 0.0);
            }
        }
    }
    // Grab inact-act and act-vir orbs
    // C is actually (nmo x nmo)
    std::shared_ptr<Matrix> Cia = C->get_block({dim_zero, nmopi_}, {dim_zero, noccpi});
    std::shared_ptr<Matrix> Cav = C->get_block({dim_zero, nmopi_}, {doccpi_, doccpi_ + virpi});

    // Back transform MOgradient
    std::shared_ptr<Matrix> gradient = linalg::triplet(Cia, MOgradient, Cav, false, false, true);

    return gradient;
}



} // End namespaces
