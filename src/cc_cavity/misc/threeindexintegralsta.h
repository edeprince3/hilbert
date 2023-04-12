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

#ifndef THREEINDEXINTEGRALSTA_H
#define THREEINDEXINTEGRALSTA_H
#include <psi4/libpsi4util/process.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/basisset.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libmints/sieve.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/psifiles.h>
namespace hilbert {

/// Grab three-index integrals from CD basis
double* ThreeIndexIntegrals(std::shared_ptr<psi::Wavefunction> ref, size_t &nQ, long int memory);

}

#endif
