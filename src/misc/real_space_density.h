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

#ifndef REAL_SPACE_DENSITY_H
#define REAL_SPACE_DENSITY_H

#define PSIF_DCC_QMO          268
#define PSIF_V2RDM_CHECKPOINT 269
#define PSIF_V2RDM_D2AA       270
#define PSIF_V2RDM_D2AB       271
#define PSIF_V2RDM_D2BB       272
#define PSIF_V2RDM_D3AAA      273
#define PSIF_V2RDM_D3AAB      274
#define PSIF_V2RDM_D3BBA      275
#define PSIF_V2RDM_D3BBB      276
#define PSIF_V2RDM_D1A        277
#define PSIF_V2RDM_D1B        278

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string>

#include "psi4/libmints/wavefunction.h"

// for reading integrals from disk
#include <psi4/libiwl/iwl.h>

// for dft
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"

// for grid
#include "psi4/libfock/points.h"
#include "psi4/libfock/cubature.h"

#include "psi4/psi4-dec.h"
#include <psi4/psifiles.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libpsi4util/PsiOutStream.h>

// tpdm and opdm structs live here
#include <v2rdm_casscf/v2rdm_solver.h>

using namespace psi;

namespace hilbert{ 

class RealSpaceDensity: public Wavefunction{

  public:

    RealSpaceDensity(std::shared_ptr<psi::Wavefunction> reference_wavefunction,Options & options);
    ~RealSpaceDensity();

    /// return grid points (x)
    std::shared_ptr<Vector> grid_x() { return grid_x_; }

    /// return grid points (y)
    std::shared_ptr<Vector> grid_y() { return grid_y_; }

    /// return grid points (z)
    std::shared_ptr<Vector> grid_z() { return grid_z_; }

    /// return grid weights (w)
    std::shared_ptr<Vector> grid_w() { return grid_w_; }

    /// return xc hole on grid 
    std::shared_ptr<Vector> xc_hole(double x, double y, double z);

    /// return on-top pair density (pi) on grid 
    std::shared_ptr<Vector> pi() { 
        BuildPiFromDisk();
        return pi_; 
    }

    /// return density (rho_a + rho_b) on grid 
    std::shared_ptr<Vector> rho() { return rho_; }

    /// return density (rho_a) on grid 
    std::shared_ptr<Vector> rho_a() { return rho_a_; }

    /// return density (rho_b) on grid 
    std::shared_ptr<Vector> rho_b() { return rho_b_; }

    /// return derivative of alpha-spin density with respect to x on grid 
    std::shared_ptr<Vector> rho_a_x() { return rho_a_x_; }

    /// return derivative of alpha-spin density with respect to y on grid 
    std::shared_ptr<Vector> rho_a_y() { return rho_a_y_; }

    /// return derivative of alpha-spin density with respect to z on grid 
    std::shared_ptr<Vector> rho_a_z() { return rho_a_z_; }

    /// return derivative of beta-spin density with respect to x on grid 
    std::shared_ptr<Vector> rho_b_x() { return rho_b_x_; }

    /// return derivative of beta-spin density with respect to y on grid 
    std::shared_ptr<Vector> rho_b_y() { return rho_b_y_; }

    /// return derivative of beta-spin density with respect to z on grid 
    std::shared_ptr<Vector> rho_b_z() { return rho_b_z_; }

    /// return the alpha opdm
    std::shared_ptr<Matrix> Da() { return Da_; }

    /// return the beta opdm
    std::shared_ptr<Matrix> Db() { return Db_; }

    void common_init();

    /// build real-space density (rho) on grid, from 1RDM on disk
    void BuildRhoFromDisk();

    virtual bool same_a_b_orbs() const { return same_a_b_orbs_; }
    virtual bool same_a_b_dens() const { return same_a_b_dens_; }

  protected:

    /// nonzero elements of alpha opdm
    std::vector<opdm> opdm_a_;

    /// nonzero elements of beta opdm
    std::vector<opdm> opdm_b_;

    /// dft potential object
    std::shared_ptr<VBase> potential_;

    /// points function
    std::shared_ptr<PointFunctions> points_func_;

    /// maximum number of functions in a block of the phi matrix
    int max_functions_;

    /// maximum number of grid points in a block of the phi matrix
    int max_points_;

    /// orbital symmetry
    int * symmetry_;

    /// geminals, by symmetry
    std::vector < std::vector < std::pair<int,int> > > gems_;

    /// reference energy
    double reference_energy_;

    /// offset for orbitals in each irrep
    int *pitzer_offset_;

    /// available memory
    long int available_memory_;

    /// number of grid points_
    long int phi_points_;

    /// phi matrix (MO)
    std::shared_ptr<Matrix> super_phi_;

    /// d phi / dx matrix (MO)
    std::shared_ptr<Matrix> super_phi_x_;

    /// d phi / dy matrix (MO)
    std::shared_ptr<Matrix> super_phi_y_;

    /// d phi / dz matrix (MO)
    std::shared_ptr<Matrix> super_phi_z_;

    /// grid x values
    std::shared_ptr<Vector> grid_x_;

    /// grid y values
    std::shared_ptr<Vector> grid_y_;

    /// grid z values
    std::shared_ptr<Vector> grid_z_;

    /// grid weights
    std::shared_ptr<Vector> grid_w_;
    
    /// get grid coordinates and weights
    void GetGridInfo();

    /// a function to build phi/phi_x/... note: the orbital labels are in the AO basis (no symmetry)
    void BuildPhiMatrixAO(std::string phi_type, std::shared_ptr<Matrix> myphi);

    /// transform the orbital labels in phi/phi_x/... from the AO to the MO basis
    void TransformPhiMatrixAOMO(std::shared_ptr<Matrix> phi_in, std::shared_ptr<Matrix> phi_out);

    /// read 2-RDM from disk and build on-top pair density
    void BuildPiFromDisk();

    /// exchange-correlation hole
    std::shared_ptr<Vector> xc_hole_;

    /// alpha-spin density
    std::shared_ptr<Vector> rho_a_;

    /// beta-spin density
    std::shared_ptr<Vector> rho_b_;

    /// total density
    std::shared_ptr<Vector> rho_;
 
    /// alpha-spin density gradient x
    std::shared_ptr<Vector> rho_a_x_;

    /// beta-spin density gradient x
    std::shared_ptr<Vector> rho_b_x_;

    /// alpha-spin density gradient y
    std::shared_ptr<Vector> rho_a_y_;

    /// beta-spin density gradient y
    std::shared_ptr<Vector> rho_b_y_;

    /// alpha-spin density gradient z
    std::shared_ptr<Vector> rho_a_z_;

    /// beta-spin density gradient z
    std::shared_ptr<Vector> rho_b_z_;
    
    /// the on-top pair density
    std::shared_ptr<Vector> pi_;

    /// x-component of the gradient of the on-top pair density
    std::shared_ptr<Vector> pi_x_;

    /// y-component of the gradient of the on-top pair density
    std::shared_ptr<Vector> pi_y_;

    /// z-component of the gradient of the on-top pair density
    std::shared_ptr<Vector> pi_z_;

    /// build spin densities and gradients using only non-zero elements of OPDM
    void BuildRhoFast();

    /// build on-top pair density using only non-zero elements of TPDM
    void BuildPiFast(std::vector<tpdm> D2ab, int nab);

    /// build exchange correlation hole
    void BuildExchangeCorrelationHole(size_t p);

};

} // end of namespaces

#endif
