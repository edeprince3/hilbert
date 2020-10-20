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

    /// alpha + beta Fock matrix
    std::shared_ptr<Matrix> F_;

    /// alpha + beta dipole x
    std::shared_ptr<Matrix> Dipole_x_;

    /// alpha + beta dipole y
    std::shared_ptr<Matrix> Dipole_y_;

    /// alpha + beta dipole z
    std::shared_ptr<Matrix> Dipole_z_;

    /// alpha + beta core hamiltonian matrix
    std::shared_ptr<Matrix> H_;

    /// alpha + beta extra one-electron terms introduced by cavity
    std::shared_ptr<Matrix> oe_cavity_terms_;

    /// number of auxiliary basis functions
    size_t nQ_;

    /// temporary container 1
    double * tmp1_;

    /// temporary container 2
    double * tmp2_;

    /// temporary container 3
    double * tmp3_;

    /// the number of cc amplitudes
    size_t ccamps_dim_;

    /// cc amplitudes (all)
    double * ccamps_;

    /// t amplitudes (doubles)
    double * t2_;

    /// t amplitudes (singles)
    double * t1_;

    /// do account for u0?
    bool include_u0_;

    /// do account for u1?
    bool include_u1_;

    /// do account for u2?
    bool include_u2_;

    /// u amplitudes (0)
    double * u0_;

    /// u amplitudes (singles)
    double * u1_;

    /// u amplitudes (doubles)
    double * u2_;

    /// residual (all)
    double * residual_;

    /// residual (t doubles)
    double * rt2_;

    /// residual (t singles)
    double * rt1_;

    /// residual (u doubles)
    double * ru2_;

    /// residual (u singles)
    double * ru1_;

    /// residual (u 0)
    double * ru0_;

    /// orbital energies (full a+b list)
    double * epsilon_;

    /// the DIIS solver
    std::shared_ptr<DIIS> diis;

    /// solve CC equations
    double cc_iterations();

    /// build total residual
    void residual();

    /// build t1 part of residual
    void residual_t1();

    /// build t2 part of residual
    void residual_t2();

    /// build u1 part of residual
    void residual_u1();

    /// build u2 part of residual
    void residual_u2();

    /// build u0 part of the residual
    void residual_u0();

    /// update t amplitudes
    double update_amplitudes();

    /// evaluate correlation energy
    double correlation_energy();

    /// build t1-transformed integrals
    double t1_transformation();

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

    /// <ab||ij>
    double * eri_abij_;

    /// <ia||jb>
    double * eri_iajb_;

    /// <jk||ia>
    double * eri_jkia_;

    /// <ia||jk>
    double * eri_iajk_;

    /// <ai||bc>
    double * eri_aibc_;

    /// <ab||ic>
    double * eri_abic_;

};

} // End namespaces

#endif 
