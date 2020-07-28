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

#ifndef JELLIUM_SCF_H
#define JELLIUM_SCF_H

#include <vector>

#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/matrix.h>

#include "jellium_integrals.h"

using namespace psi;

namespace hilbert{

struct cis_transition {
    int i;
    int a;
    int hi;
    int ha;
};

class Jellium_SCFSolver{

  public:

    Jellium_SCFSolver(Options & options);
    ~Jellium_SCFSolver();
    double compute_energy();

    /// construct CIS sigma vectors for out-of-core davidson solver
    void CIS_evaluate_sigma(size_t N, size_t maxdim, double ** bmat, double ** sigma);

  protected:

    void CIS_direct();

    void CIS_slow();

    /// Options object
    Options& options_;    

    /// finite jellium integrals
    std::shared_ptr<JelliumIntegrals> jelly_;

    /// number of orbitals
    int nso_;

    /// number of irreps
    int nirrep_;

    /// number of orbitals per irrep
    int * nsopi_;

    /// total number of electrons
    int nelectron_;

    /// number of doubly occupied orbitals per irrep
    int * doccpi_;

    /// length of box
    double boxlength_;

    /// scaling factor for integrals based on box length
    double Lfac_;

    /// kinetic energy integrals
    std::shared_ptr<Matrix> T_;

    /// potential energy integrals
    std::shared_ptr<Matrix> V_;

    /// so/mo transformation matrix
    std::shared_ptr<Matrix> Ca_;

    /// so-basis density matrix
    std::shared_ptr<Matrix> Da_;

    /// so-basis fock matrix
    std::shared_ptr<Matrix> Fa_;

    /// build coulomb matrix
    void build_J(std::shared_ptr<Matrix> Da, std::shared_ptr<Matrix> Ja);

    /// build exchange matrix
    void build_K(std::shared_ptr<Matrix> Da, std::shared_ptr<Matrix> Ka);

    /// list of CIS excitations
    std::vector<cis_transition> cis_transition_list_;

    /// update occupations
    void update_occupations(std::shared_ptr<Vector> epsilon_a);

    /// print occupations
    void print_occupations();

};

} // End namespaces



#endif
