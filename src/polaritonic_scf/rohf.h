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

#ifndef POLARITONIC_ROHF_H
#define POLARITONIC_ROHF_H

#include "hf.h"

using namespace psi;

namespace hilbert{ 

class PolaritonicROHF: public PolaritonicHF {

  public:

    PolaritonicROHF(std::shared_ptr<Wavefunction> reference_wavefunction,Options & options);

    ~PolaritonicROHF();

    void common_init();

    double compute_energy();

  protected:

    /// form effective fock matrix
    void form_Feff(std::shared_ptr<Matrix> Feff);

    /// form ROHF orbital gradient
    std::shared_ptr<Matrix> OrbitalGradientROHF(std::shared_ptr<Matrix> F, std::shared_ptr<Matrix> C);


};

} // End namespaces

#endif 
