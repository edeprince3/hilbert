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

#include <string>

#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsi4util/PsiOutStream.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/molecule.h>

#include "orbital_optimizer.h"

using namespace psi;

namespace hilbert{

OrbitalOptimizer::OrbitalOptimizer(std::shared_ptr<Wavefunction> reference_wavefunction, Options & options){

    is_energy_converged_ = false;
    is_gradient_converged_ = false;

    nirrep_ = reference_wavefunction->nirrep();
    nmo_    = reference_wavefunction->nmo();
    doccpi_ = reference_wavefunction->doccpi();
    frzcpi_ = reference_wavefunction->frzcpi();
    frzvpi_ = reference_wavefunction->frzvpi();
    nmopi_  = reference_wavefunction->nmopi();

    // restricted doubly occupied orbitals per irrep (optimized)
    rstcpi_   = (int*)malloc(nirrep_*sizeof(int));
    memset((void*)rstcpi_,'\0',nirrep_*sizeof(int));

    // restricted unoccupied occupied orbitals per irrep (optimized)
    rstvpi_   = (int*)malloc(nirrep_*sizeof(int));
    memset((void*)rstvpi_,'\0',nirrep_*sizeof(int));

    // active orbitals per irrep:
    amopi_    = (int*)malloc(nirrep_*sizeof(int));
    memset((void*)amopi_,'\0',nirrep_*sizeof(int));

    if (options["FROZEN_DOCC"].has_changed()) {
        if (options["FROZEN_DOCC"].size() != nirrep_) {
            throw PsiException("The FROZEN_DOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzcpi_[h] = options["FROZEN_DOCC"][h].to_double();
        }
    }
    if (options["RESTRICTED_DOCC"].has_changed()) {
        if (options["RESTRICTED_DOCC"].size() != nirrep_) {
            throw PsiException("The RESTRICTED_DOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            rstcpi_[h] = options["RESTRICTED_DOCC"][h].to_double();
        }
    }
    if (options["RESTRICTED_UOCC"].has_changed()) {
        if (options["RESTRICTED_UOCC"].size() != nirrep_) {
            throw PsiException("The RESTRICTED_UOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            rstvpi_[h] = options["RESTRICTED_UOCC"][h].to_double();
        }
    }
    if (options["FROZEN_UOCC"].has_changed()) {
        if (options["FROZEN_UOCC"].size() != nirrep_) {
            throw PsiException("The FROZEN_UOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzvpi_[h] = options["FROZEN_UOCC"][h].to_double();
        }
    }
    // user could specify active space with ACTIVE array
    if ( options["ACTIVE"].has_changed() ) {
        if (options["ACTIVE"].size() != nirrep_) {
            throw PsiException("The ACTIVE array has the wrong dimensions_",__FILE__,__LINE__);
        }

        // warn user that active array takes precedence over restricted_uocc array
        if (options["RESTRICTED_UOCC"].has_changed()) {
            outfile->Printf("\n");
            outfile->Printf("    <<< WARNING!! >>>\n");
            outfile->Printf("\n");
            outfile->Printf("    The ACTIVE array takes precedence over the RESTRICTED_UOCC array.\n");
            outfile->Printf("    Check below whether your active space was correctly specified.\n");
            outfile->Printf("\n");
        }

        // overwrite rstvpi_ array using the information in the frozen_docc,
        // restricted_docc, active, and frozen_virtual arrays.  start with nso total
        // orbitals and let the linear dependency check below adjust the spaces as needed
        for (int h = 0; h < nirrep_; h++) {
            amopi_[h]  = options["ACTIVE"][h].to_double();
            rstvpi_[h] = nmopi_[h] - frzcpi_[h] - rstcpi_[h] - frzvpi_[h] - amopi_[h];
        }
    }

    // now, total number of orbitals in each class:
    nfrzc_ = 0;
    nrstc_ = 0;
    amo_   = 0;
    nrstv_ = 0;
    nfrzv_ = 0;
    for (int h = 0; h < nirrep_; h++) {
        nfrzc_   += frzcpi_[h];
        nrstc_   += rstcpi_[h];
        nrstv_   += rstvpi_[h];
        nfrzv_   += frzvpi_[h];
        amo_     += nmopi_[h]-frzcpi_[h]-rstcpi_[h]-rstvpi_[h]-frzvpi_[h];
        amopi_[h] = nmopi_[h]-frzcpi_[h]-rstcpi_[h]-rstvpi_[h]-frzvpi_[h];
    }

    // lastly, determine orbital symmetries (in energy order)
    orbital_symmetries_ = (int*)malloc(nmo_*sizeof(int));
    memset((void*)orbital_symmetries_,'\0',nmo_*sizeof(int));

    bool * skip = (bool*)malloc(nmo_*sizeof(bool));
    memset((void*)skip,'\0',nmo_*sizeof(bool));

    // frozen core
    for (int i = 0; i < nfrzc_; i++){
        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) { 
            for (int j = 0; j < frzcpi_[h]; j++){
                if ( reference_wavefunction->epsilon_a()->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = reference_wavefunction->epsilon_a()->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            me += nmopi_[h] - frzcpi_[h] - frzvpi_[h];
        }
        if ( imin < 0 ) {
            throw PsiException("cannot determine energy-to-pitzer order (frzc)",__FILE__,__LINE__);
        }
        skip[imin] = true;
        orbital_symmetries_[i] = isym;
    }

    // restricted doubly occupied
    for (int i = nfrzc_; i < nrstc_ + nfrzc_; i++){

        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) {
            me += frzcpi_[h];
            for (int j = frzcpi_[h]; j < rstcpi_[h] + frzcpi_[h]; j++){
                if ( reference_wavefunction->epsilon_a()->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = reference_wavefunction->epsilon_a()->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            me += amopi_[h] + rstvpi_[h];// + frzvpi_[h];
        }
        if ( imin < 0 ) {
            throw PsiException("cannot determine energy-to-pitzer order (rstc)",__FILE__,__LINE__);
        }
        skip[imin] = true;
        orbital_symmetries_[i] = isym;
    }
    // active
    for (int i = nrstc_ + nfrzc_; i < amo_ + nrstc_ + nfrzc_; i++){

        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) {
            me += rstcpi_[h] + frzcpi_[h];
            for (int j = rstcpi_[h] + frzcpi_[h]; j < rstcpi_[h] + frzcpi_[h]+amopi_[h]; j++){
                if ( reference_wavefunction->epsilon_a()->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = reference_wavefunction->epsilon_a()->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            me += rstvpi_[h];// + frzvpi_[h];
        }
        if ( imin < 0 ) {
            throw PsiException("cannot determine energy-to-pitzer order (active)",__FILE__,__LINE__);
        }
        skip[imin] = true;
        orbital_symmetries_[i] = isym;
    }
    // restricted virtual
    for (int i = amo_ + nrstc_ + nfrzc_; i < nmo_ - nfrzv_ ; i++){

        int me     = 0;
        double min = 1.0e99;
        int imin   = -999;
        int isym   = -999;
        for (int h = 0; h < nirrep_; h++) {
            me += rstcpi_[h] + frzcpi_[h] + amopi_[h];
            for (int j = rstcpi_[h] + frzcpi_[h] + amopi_[h]; j < nmopi_[h] - frzvpi_[h] ; j++){
                if ( reference_wavefunction->epsilon_a()->pointer(h)[j] < min ) {
                    if ( !skip[me] ) {
                        min = reference_wavefunction->epsilon_a()->pointer(h)[j];
                        imin = me;
                        isym = h;
                    }
                }
                me++;
            }
            //me += frzvpi_[h];
        }
        if ( imin < 0 ) {
            throw PsiException("cannot determine energy-to-pitzer order (rstv)",__FILE__,__LINE__);
        }
        skip[imin] = true;
        orbital_symmetries_[i] = isym;
    }

    free(skip);

    // print orbitals per irrep in each space
    outfile->Printf("  ==> Active space details <==\n");
    outfile->Printf("\n");
    //outfile->Printf("        Freeze core orbitals?                   %5s\n",nfrzc_ > 0 ? "yes" : "no");
    outfile->Printf("        Number of frozen core orbitals:         %5i\n",nfrzc_);
    outfile->Printf("        Number of restricted occupied orbitals: %5i\n",nrstc_);
    outfile->Printf("        Number of active orbitals:              %5i\n",amo_);
    outfile->Printf("        Number of restricted virtual orbitals:  %5i\n",nrstv_);
    outfile->Printf("        Number of frozen virtual orbitals:      %5i\n",nfrzv_);
    outfile->Printf("\n");

    std::vector<std::string> labels = reference_wavefunction->molecule()->irrep_labels();
    outfile->Printf("        Irrep:           ");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4s",labels[h].c_str());
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" \n");
    outfile->Printf(" \n");

    outfile->Printf("        frozen_docc     [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",frzcpi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("        restricted_docc [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",rstcpi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("        active          [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",amopi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("        restricted_uocc [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",rstvpi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("        frozen_uocc     [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",frzvpi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("\n");

    // options
    // TODO: some options are missing:
    //    ORBOPT_ONE_STEP (should be handled outside of this class)
    //    ORBOPT_FREQUENCY (should be handled outside of this class)
    //    ORBOPT_NUM_DIIS_VECTORS (is diis implemented?)

    e_convergence_           = options.get_double("ORBOPT_ENERGY_CONVERGENCE");
    g_convergence_           = options.get_double("ORBOPT_GRADIENT_CONVERGENCE");
    active_active_rotations_ = options.get_bool("ORBOPT_ACTIVE_ACTIVE_ROTATIONS");
    maxiter_                 = options.get_int("ORBOPT_MAXITER");
    write_                   = options.get_bool("ORBOPT_WRITE");
    exact_diagonal_hessian_  = options.get_bool("ORBOPT_EXACT_DIAGONAL_HESSIAN");
    algorithm_               = options.get_str("ORBOPT_ALGORITHM");

    outfile->Printf("  ==> Orbital optimization parameters <==\n");
    outfile->Printf("\n");
    outfile->Printf("        g_convergence:               %12.3le\n",g_convergence_);
    outfile->Printf("        e_convergence:               %12.3le\n",e_convergence_);
    outfile->Printf("        maximum iterations:          %12i\n",maxiter_);
    outfile->Printf("        exact diagonal Hessian:      %12s\n",exact_diagonal_hessian_ ? "true" : "false");
    outfile->Printf("        print iteration info:        %12s\n",write_ ? "true" : "false");
    outfile->Printf("        algorithm:                   %12s\n",algorithm_.c_str());
    outfile->Printf("        active-active rotations:     %12s\n",active_active_rotations_ ? "true" : "false");
    outfile->Printf("\n");

}

OrbitalOptimizer::~OrbitalOptimizer(){

    free(rstcpi_);
    free(rstvpi_);
    free(amopi_);
    free(orbital_symmetries_);

}

void OrbitalOptimizer::optimize_orbitals(double * d2, double * d1, double * tei, double * oei, double * transformation_matrix){
    throw PsiException("optimize_orbitals not yet implemented",__FILE__,__LINE__);
}

void OrbitalOptimizer::get_lagrangian(std::shared_ptr<Matrix> Lagrangian){
    throw PsiException("get_lagrangian not yet implemented",__FILE__,__LINE__);
}

void OrbitalOptimizer::get_hessian(std::shared_ptr<Matrix> Hessian){
    throw PsiException("get_hessian not yet implemented",__FILE__,__LINE__);
}



}

