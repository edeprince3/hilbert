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

#ifndef POLARITONIC_RDFRPA_H
#define POLARITONIC_RDFRPA_H

#include "hf.h"
#include "rtddft.h"
#include <psi4/lib3index/dfhelper.h>

using namespace psi;

namespace hilbert{ 

class PolaritonicRDFRPA: public PolaritonicHF {

  public:

    PolaritonicRDFRPA(std::shared_ptr<Wavefunction> reference_wavefunction,Options & options,std::shared_ptr<Wavefunction> dummy_wfn);

    ~PolaritonicRDFRPA();

    void common_init(std::shared_ptr<Wavefunction> dummy_wfn);
    // batched version
    void common_init_batched(std::shared_ptr<Wavefunction> dummy_wfn);

    // dispatch function for energy computation
    double compute_energy();
    // unbatched version
    double compute_energy_unbatched();
    // batched version
    double compute_energy_batched();

    // batchsize
    int determine_optimal_batch_size();

  private:

    std::vector<double> project_dipole_into_aux(const std::vector<double>& mu_ai);

  protected:

    std::shared_ptr<Vector> build_diag_matrix();

    std::shared_ptr<Matrix> build_grid();

    std::shared_ptr<Matrix> build_minimax_grid(double tol);

    std::shared_ptr<Matrix> tst_;
    
    std::shared_ptr<psi::DFHelper> df_helper_;

    long int nQ_;
    long int nQ_eff_;
    long int ia_batch_size_;
    long int o_; 
    long int v_;

    // control output (for debugging, to be done properly)
    long int n_output_ = 2;

    // number of grid points for quadrature
    long int n_grid_points_ = 60;

    // scaling factor for quadrature
    double scal_fac_ = 1.0;

    // difference of orbital energies
    double * diag_;
    // ov block of the 3 index integrals
    double * siap_;
    // temporary container for the 3 index integrals
    double * tmp1_;
    // matrix Q accumulated over quadrature points
    double * Qmat_;
    // integer array for dsyevd
    int * int1_;
};

} // End namespaces

#endif 
