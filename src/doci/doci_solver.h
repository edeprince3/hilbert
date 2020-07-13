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

#ifndef DOCI_SOLVER_H
#define DOCI_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <bitset>

#include <psi4/libiwl/iwl.h>
#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libqt/qt.h>

#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>

#include <misc/hilbert_psifiles.h>

using namespace psi;

namespace hilbert{

class DOCISolver: public Wavefunction{
  public:
    DOCISolver(SharedWavefunction reference_wavefunction,Options & options);
    ~DOCISolver();
    double compute_energy();

    // construct sigma vectors for out-of-core davidson solver
    void BuildSigma(size_t N, size_t maxdim, double ** b, double ** sigma);

    // construct sigma vectors for out-of-core davidson solver
    void BuildSigmaFast(size_t N, size_t maxdim, double ** b, double ** sigma);

    // construct a given element of the hamiltonian
    double HamiltonianElement(size_t i, size_t j);


  protected:

    virtual bool same_a_b_orbs() const { return same_a_b_orbs_; }
    virtual bool same_a_b_dens() const { return same_a_b_dens_; }

    void common_init();

    /// bit strings representing configurations
    std::bitset<64> * configurations_;

    /// which configurations are connected to which?
    size_t ** configuration_map_;

    /// what orbitals differ between connected configurations?
    size_t *** configuration_map_indices_;

    /// 1-RDM
    double * d1_;

    /// 2-RDM (ab)
    double * d2ab_;

    /// 2-RDM (aa)
    double * d2aa_;

    /// build and diagonalize Hamiltonian
    double DiagonalizeHamiltonian(double * ci_wfn,size_t & iter);

    /// diagonals of Hamiltonian
    std::shared_ptr<Vector> Hdiag_;

    /// construct 1- and 2-RDM
    void BuildRDMs(double * ci_wfn,bool print);

    /// optimize orbitals
    double RotateOrbitals();

    // total number of configurations
    size_t n_;

    /// generate all doubly occupied configurations
    size_t GenerateConfigurations(size_t val);

    /// are we using density fitting?
    bool is_df_;

    /// convergence in energy
    double e_convergence_;

    /// convergence in amplitudes
    double r_convergence_;

    /// three-index integrals in the MO basis
    double * Qmo_;

    /// number of auxiliary basis functions
    long int nQ_;

    /// SCF energy
    double escf_;

    /// nuclear repulsion energy
    double enuc_;

    /// full space of integrals for MO gradient / Hessian, blocked by symmetry
    double * tei_;
    long int tei_dim_;

    /// (ii|jj)
    double * eint1_;

    /// (ij|ji)
    double * eint2_;

    double * oei_;
    int oei_dim_;
    int d2_dim_;
    double * d2_;

    /// pack active-space spin-blocked density into spatial density
    void PackSpatialDensity();

    /// orbital optimization stuff:

    double * orbopt_transformation_matrix_;
    double * orbopt_data_;
    char * orbopt_outfile_;
    bool orbopt_converged_;

    /// wall time for orbital optimization
    double orbopt_time_;

    /// maximum number of macroiterations
    size_t maxiter_;

    /// write full 2RDM to disk
    void WriteTPDM();

    /// mo-mo transformation matrix
    SharedMatrix newMO_;

    /// print 1- and 2-RDM elements
    void print_rdms();

};

}
#endif

