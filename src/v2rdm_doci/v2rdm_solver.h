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

#ifndef V2RDM_SOLVER_H
#define V2RDM_SOLVER_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

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

namespace psi{ namespace v2rdm_doci{

class v2RDMSolver: public Wavefunction{
  public:
    v2RDMSolver(SharedWavefunction reference_wavefunction,Options & options);
    ~v2RDMSolver();
    void common_init();
    double compute_energy();
    virtual bool same_a_b_orbs() const { return same_a_b_orbs_; }
    virtual bool same_a_b_dens() const { return same_a_b_dens_; }

    // public methods
    void cg_Ax(long int n,SharedVector A, SharedVector u);

  protected:

    /// constrain T1 to be positive semidefinite?
    bool constrain_t1_;

    /// constrain T2 to be positive semidefinite?
    bool constrain_t2_;

    /// constrain Q2 to be positive semidefinite?
    bool constrain_q2_;

    /// constrain G2 to be positive semidefinite?
    bool constrain_g2_;

    /// constrain D3 to be positive semidefinite?
    bool constrain_d3_;

    /// symmetry product table:
    int * table;

    /// returns symmetry product for two orbitals
    int SymmetryPair(int i, int j);

    /// returns symmetry product for four orbitals
    int TotalSym(int i, int j,int k, int l);

    int * symmetry;
    int * symmetry_full;
    int * symmetry_really_full;
    int * energy_to_pitzer_order;
    int * energy_to_pitzer_order_really_full;
    int * symmetry_energy_order;
    int * pitzer_offset;
    int * pitzer_offset_full;

    /// active-space geminals for each symmetry:
    std::vector < std::vector < std::pair<int,int> > > gems;

    /// total number of active molecular orbitals
    int amo_;

    /// total number of frozen core orbitals
    int nfrzc_;

    /// total number of frozen virtual orbitals
    int nfrzv_;

    /// total number of restricted doubly occupied orbitals
    int nrstc_;

    /// total number of restricted unoccupied orbitals
    int nrstv_;

    /// active molecular orbitals per irrep
    int * amopi_;

    /// restricted core orbitals per irrep.  these will be optimized optimized
    int * rstcpi_;

    /// restricted virtual orbitals per irrep.  these will be optimized optimized
    int * rstvpi_;

    /// total number of constraints (dimension of dual solution vector)
    long int nconstraints_;

    /// total number of variables (dimension of primal solution vector)
    long int dimx_;

    /// number of auxilliary basis functions
    long int nQ_;

    /// read three-index integrals and transform them to MO basis
    void ThreeIndexIntegrals();

    /// three-index integral buffer
    double * Qmo_;

    /// grab one-electron integrals (T+V) in MO basis
    SharedMatrix GetOEI();

    /// offsets
    int d1off_;
    int d2s2off_;
    int d2s0off_;
    int q2s2off_;
    int q2s0off_;
    int g2s2off_;
    int g2s0off_;
    int t1s3aaaoff_;
    int t1s3aaboff_;
    int t1s1off_;
    int t2s3off_;
    int t2s1off_;
    int t2off_;
    int d3s3off_;
    int d3s1off_;

    /// convergence in primal energy
    double e_convergence_;

    /// convergence in primal and dual error
    double r_convergence_;

    /// convergence in conjugate gradient solver
    double cg_convergence_;

    /// maximum number of boundary-point SDP (outer) iterations
    int maxiter_;

    /// maximum number of conjugate gradient (inner) iterations
    int cg_maxiter_;

    /// standard vector of dimensions of each block of primal solution vector
    std::vector<int> dimensions_;

    int offset;

    // mapping arrays with abelian symmetry
    void BuildBasis();
    int * full_basis;

    /// mapping arrays with symmetry:
    int * gems_ab;
    int * gems_aa;
    int * gems_00;
    int * gems_full;
    int * gems_plus_core;
    int *** bas_ab_sym;
    int *** bas_aa_sym;
    int *** bas_00_sym;
    int *** bas_full_sym;
    int *** bas_really_full_sym;
    int *** ibas_ab_sym;
    int *** ibas_aa_sym;
    int *** ibas_00_sym;
    int *** ibas_full_sym;
    int *** ibas_really_full_sym;

    void PrintHeader();

    /// grab one- and two-electron integrals
    void GetIntegrals();

    /// read two-electron integrals from disk
    void GetTEIFromDisk();

    // read teis from disk:
    void ReadAllIntegrals(iwlbuf *Buf);

    /// grab a specific two-electron integral
    double TEI(int i, int j, int k, int l, int h);

    void BuildConstraints();

    void Guess();

    void bpsdp_Au(SharedVector A, SharedVector u);
    void D2_constraints_Au(SharedVector A,SharedVector u);
    void Q2_constraints_Au(SharedVector A,SharedVector u);
    void G2_constraints_Au(SharedVector A,SharedVector u);
    void D3_constraints_Au(SharedVector A,SharedVector u);
    void T1_constraints_Au(SharedVector A,SharedVector u);
    void T2_constraints_Au(SharedVector A,SharedVector u);

    void bpsdp_ATu(SharedVector A, SharedVector u);
    void D2_constraints_ATu(SharedVector A,SharedVector u);
    void Q2_constraints_ATu(SharedVector A,SharedVector u);
    void G2_constraints_ATu(SharedVector A,SharedVector u);
    void D3_constraints_ATu(SharedVector A,SharedVector u);
    void T1_constraints_ATu(SharedVector A,SharedVector u);
    void T2_constraints_ATu(SharedVector A,SharedVector u);

    /// SCF energy
    double escf_;

    /// nuclear repulsion energy
    double enuc_;

    double tau, mu, ed, ep;

    //vectors
    SharedVector Ax;     // vector to hold A . x
    SharedVector ATy;    // vector to hold A^T . y
    SharedVector c;      // 1ei and 2ei of bpsdp
    SharedVector y;      // dual solution
    SharedVector b;      // constraint vector
    SharedVector x;      // primal solution
    SharedVector z;      // second dual solution

    void Update_xz();
    void Update_xz_nonsymmetric();

    void NaturalOrbitals();
    void MullikenPopulations();

    // multiplicity
    int multiplicity_;

    /// full space of integrals for MO gradient / Hessian, blocked by symmetry
    double * tei_full_sym_;
    long int tei_full_dim_;

    double * oei_full_sym_;
    int oei_full_dim_;

    /// full space D2, blocked by symmetry
    double * d2_plus_core_sym_;
    int d2_plus_core_dim_;

    /// active space spatial D2, blocked by symmetry
    double * d2_act_spatial_sym_;
    int d2_act_spatial_dim_;

    /// active space spatial D1, blocked by symmetry
    double * d1_act_spatial_sym_;
    int d1_act_spatial_dim_;

    /// pack active-space spin-blocked density into spatial density
    void PackSpatialDensity();

    /// repack rotated full-space integrals into active-space integrals
    void RepackIntegrals();

    /// compute frozen core energy and adjust oeis
    void FrozenCoreEnergy();

    /// function to rotate orbitals
    void RotateOrbitals();

    double * orbopt_transformation_matrix_;
    double * orbopt_data_;
    char * orbopt_outfile_;
    bool orbopt_converged_;

    /// are we using 3-index integrals?
    bool is_df_;

    /// write primal, dual, and orbitals to a checkpoint file
    void WriteCheckpointFile();

    /// read primal, dual, and orbitals from a checkpoint file
    void ReadFromCheckpointFile();

    /// read orbitals from a checkpoint file
    void ReadOrbitalsFromCheckpointFile();

    /// wall time for microiterations
    double iiter_time_;

    /// wall time for macroiterations
    double oiter_time_;

    /// wall time for orbital optimization
    double orbopt_time_;

    /// total number of microiterations
    long int iiter_total_;

    /// total number of macroiterations
    long int oiter_total_;

    /// total number of orbital optimization
    long int orbopt_iter_total_;

    /// write full 2RDM to disk
    void WriteTPDM();

    /// write full 2RDM to disk in IWL format
    void WriteTPDM_IWL();

    /// write active-active-active-active 2RDM to disk
    //void WriteActiveTPDM();

    /// read 2RDM from disk
    void ReadTPDM();

    /// write molden file
    void WriteMoldenFile();

    /// orbital lagrangian
    double * X_;

    void OrbitalLagrangian();
    void DualD1Q1();
    void PrintDual();

    /// memory available beyond what is allocated for v2RDM-DOCI
    long int available_memory_;

    /// update primal solution after semicanonicalization
    void UpdatePrimal();

    /// transform a four-index quantity from one basis to another
    void TransformFourIndex(double * inout, double * tmp, SharedMatrix trans);

    /// update ao/mo transformation matrix after orbital optimization
    void UpdateTransformationMatrix();

    /// mo-mo transformation matrix
    SharedMatrix newMO_;
};

}}
#endif

