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

#ifndef POLARITONIC_UCCSD_H
#define POLARITONIC_UCCSD_H

#include "hf.h"

#include <misc/diis.h>

using namespace psi;

namespace hilbert{ 

class PolaritonicUCCSD: public PolaritonicHF {

  public:

    PolaritonicUCCSD(std::shared_ptr<Wavefunction> reference_wavefunction,Options & options);

    ~PolaritonicUCCSD();

    void common_init();

    double compute_energy();

  protected:

    /// alpha + beta MO transformation matrix
    std::shared_ptr<Matrix> C_;

    /// number of auxiliary basis functions
    int nQ_;

    /// temporary container 1
    double * tmp1_;

    /// temporary container 2
    double * tmp2_;

    /// temporary container 3
    double * tmp3_;

    /// t amplitudes (all)
    double * tamps_;

    /// t amplitudes (doubles)
    double * t2_;

    /// t amplitudes (singles)
    double * t1_;

    /// residual (all)
    double * residual_;

    /// residual (doubles)
    double * r2_;

    /// residual (singles)
    double * r1_;

    /// orbital energies (full a+b list)
    double * epsilon_;

    /// the DIIS solver
    std::shared_ptr<DIIS> diis;

    /// build residual
    void residual();

    /// update t amplitudes
    double update_amplitudes();

    /// evaluate correlation energy
    double correlation_energy();

    /// build mo-basis electron repulsion integrals
    void build_mo_eris();

    /// generate and write SO-basis three-index integrals to disk
    void write_three_index_ints();

    /// <ij||kl>
    double * eri_ijkl_;

    /// <ab||cd>
    double * eri_abcd_;

    /// <ij||ab>
    double * eri_ijab_;

    /// <ia||jb>
    double * eri_iajb_;

    /// <ia||jk>
    double * eri_iajk_;

    /// <ab||ci>
    double * eri_aibc_;

};

} // End namespaces

#endif 
