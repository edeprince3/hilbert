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
//#include "psi4/libscf_solver/hf.h"

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

    void compute_properties();

    void build_sigma_generalized(int N, int maxdim, int L, double **Q, double **sigmah, double **sigmas);

  protected:

    std::vector<std::vector<double>> compute_first_order_response(double omega);

    void compute_polarizability(std::vector<double>X, std::vector<double>Y, double omega);

    void compute_hyperpolarizability(std::vector<std::vector<double>>amps_wx, 
                                     std::vector<std::vector<double>>amps_wy,
                                     std::vector<std::vector<double>>amps_wz,
                                     std::string type, double omega);

    void build_Au_Bu(int N, int L, double *u, double *ABu);

    void build_sigma_m(int N, int L, double *x, double *y, double *m, double *sigma_m);

    void build_gm(int N, int L, double *m, double *gm);

    std::shared_ptr<VBase> potential_;

    double * build_hamiltonian_diagonals();

    double * build_overlap_diagonals();

    // jk object
    std::shared_ptr<JK> jk_;

    // is the functional hf?
    bool is_hf_;

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

    // MO-basis dipole integrals
    std::vector<std::shared_ptr<Matrix>> mua_;
    std::vector<std::shared_ptr<Matrix>> mub_;

    // MO-basis lambda-dressed dipole integrals
    std::shared_ptr<Matrix> lambda_dressed_mua_;
    std::shared_ptr<Matrix> lambda_dressed_mub_;

};

} // End namespaces

#endif 
