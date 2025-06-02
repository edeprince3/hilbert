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

#ifndef DIISTA_H
#define DIISTA_H

#include <tiledarray.h>

namespace hilbert {

class DIISTA {

public:

    DIISTA(int n);
    ~DIISTA();

    /// write Fock matrix disk
    void WriteVector(const std::vector<TA::TArrayD*>& vector);

    /// write error vector to disk
    void WriteErrorVector(const std::vector<TA::TArrayD*>& vector);

    /// perform diis extrapolation
    void Extrapolate(std::vector<TA::TArrayD*> &vector);

    /// restart diis
    void restart();

private:

    /// stores history of Fock matrix
    std::vector<std::vector<TA::TArrayD>> hist;

    /// stores history of error vector
    std::vector<std::vector<TA::TArrayD>> histerr;

    /// determine diis expansion coefficients
    void DIISCoefficients(int nvec);

    /// Compute dot product of error vectors i and j
    double diis_dot(int i, int j);

    /// Compute norm of error vector i
    double diis_norm(int i);

    /// get index for TA object
    std::string get_index(const TA::TArrayD& array);

    /// maximum number of diis vectors
    int maxdiis_;          

    /// dimension of each diis vector
    double * diisvec_;

    /// current number of diis vectors
    int diis_iter_;

    /// The error matrix shaped as maxdiis x maxdiis
    double * errmtx;

    /// diis vector to be replaced
    int replace_diis_iter_;
};

} // end of namespace

#endif
