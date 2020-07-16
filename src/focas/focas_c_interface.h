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

#ifndef FOCAS_C_INTERFACE_H
#define FOCAS_C_INTERFACE_H

using namespace psi;

namespace hilbert{

/// after orbital optimization, update Ca/Cb matrices and repack energy-order transformation matrix as pitzer order
void UpdateTransformationMatrix(std::shared_ptr<Wavefunction> ref, std::shared_ptr<Matrix> newMO,
        std::shared_ptr<Matrix> Ca, std::shared_ptr<Matrix> Cb, double * orbopt_transformation_matrix);

}

/**
 * wrappers to Greg's Fortran functions
 */

#ifndef FC_SYMBOL
#define FC_SYMBOL 2
#endif

#if   FC_SYMBOL==1
#define F77NAME(x) x
#elif FC_SYMBOL==2
#define F77NAME(x) x##_
#endif

// for greg
extern "C" {
    void F77NAME(focas_interface)(double*jacobi_transformation_matrix,
                                   double*oei_full_sym_,
                                   int &oei_full_sym_dim,
                                   double*tei_full_sym_,
                                   long int &tei_full_sym_dim,
                                   double*d1_full_sym,
                                   int &d1_full_sym_dim,
                                   double*d2_full_sym,
                                   int &d2_full_sym_dim,
                                   int *symmetry_energy_order,
//gg -- added pointer to frzcpi_ to argument list
                                   //int *frzcpi_,
                                   int &nfrzc,
                                   int &amo_,
                                   int &nfrzv,
                                   int &nirrep_,
                                   double*jacobi_data,
                                   char*jacobi_file,
                                   double*X_);
};
inline void OrbOpt(double*jacobi_transformation_matrix,
                   double*oei_full_sym_,
                   int &oei_full_sym_dim,
                   double*tei_full_sym_,
                   long int &tei_full_sym_dim,
                   double*d1_full_sym,
                   int &d1_full_sym_dim,
                   double*d2_full_sym,
                   int &d2_full_sym_dim,
                   int *symmetry_energy_order,
//gg -- added pointer to frzcpi_ to argument list
                   //int *frzcpi_,
                   int &nfrzc,
                   int &amo_,
                   int &nfrzv,
                   int &nirrep_,
                   double*jacobi_data,
//gg -- added frzcpi_ to argument list
                   char*jacobi_file,
                   double*X_){
    //F77NAME(focas_interface)(jacobi_transformation_matrix,oei_full_sym_,oei_full_sym_dim,tei_full_sym_,tei_full_sym_dim,d1_full_sym,d1_full_sym_dim,d2_full_sym,d2_full_sym_dim,symmetry_energy_order,frzcpi_,nfrzc,amo_,nfrzv,nirrep_,jacobi_data,jacobi_file);
    F77NAME(focas_interface)(jacobi_transformation_matrix,oei_full_sym_,oei_full_sym_dim,tei_full_sym_,tei_full_sym_dim,d1_full_sym,d1_full_sym_dim,d2_full_sym,d2_full_sym_dim,symmetry_energy_order,nfrzc,amo_,nfrzv,nirrep_,jacobi_data,jacobi_file,X_);
};

#endif
