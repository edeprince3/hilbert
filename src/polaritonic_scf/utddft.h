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

#ifndef POLARITONIC_UTDDFT_H
#define POLARITONIC_UTDDFT_H

#include "hf.h"

// for dft
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libscf_solver/hf.h"

// jk object
#include <psi4/libfock/jk.h>

using namespace psi;

namespace hilbert{ 

class PolaritonicUTDDFT: public PolaritonicHF {

  public:

    PolaritonicUTDDFT(std::shared_ptr<Wavefunction> reference_wavefunction,Options & options,std::shared_ptr<Wavefunction> dummy_wfn);

    ~PolaritonicUTDDFT();

    void common_init(std::shared_ptr<Wavefunction> dummy_wfn);

    double compute_energy();

    void build_sigma_generalized(int N, int maxdim, int L, double **Q, double **sigmah, double **sigmas);

  protected:

    void build_Au_Bu(int N, int L, double *u, double *Au, double *Bu);

    void build_sigma_m(int N, int L, double *x, double *y, double *m, double *sigma_m_r, double *sigma_m_l);

    void build_gm(int N, int L, double *m, double *gm);

    std::shared_ptr<VBase> potential_;

    double * build_hamiltonian_diagonals();

    double * build_overlap_diagonals();

    // jk object
    std::shared_ptr<JK> jk_;

    // hybrid?
    bool is_x_hybrid_;

    // lrc?
    bool is_x_lrc_;

    // lrc omega?
    double x_omega_;

    // needs xc?
    bool needs_xc_;

    // hybrid alpha?
    double x_alpha_;

    /// alpha + beta dipole x
    std::shared_ptr<Matrix> Dipole_x_;

    /// alpha + beta dipole y
    std::shared_ptr<Matrix> Dipole_y_;

    /// alpha + beta dipole z
    std::shared_ptr<Matrix> Dipole_z_;

    /// alpha + beta MO transformation matrix
    std::shared_ptr<Matrix> C_;

};

} // End namespaces

#endif 
