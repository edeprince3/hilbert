/*
 *@BEGIN LICENSE
 *
 * v2RDM-CASSCF, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Copyright (c) 2014, The Florida State University. All rights reserved.
 *
 *@END LICENSE
 *
 */

#ifndef V2RDM_SOLVER_H
#define V2RDM_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <psi4/libiwl/iwl.h>
#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libqt/qt.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>

#include <focas/focas_c_interface.h>
#include <misc/hilbert_psifiles.h>
#include <misc/bpsdp_solver.h>

namespace hilbert{ 

struct tpdm {
    int i;
    int j;
    int k;
    int l;
    double value;
};

struct opdm {
    int i;
    int j;
    double value;
};

// which generalized pauli constraint?
enum GeneralizedPauliConstraint {
    GeneralizedPauli_3_6,
    GeneralizedPauli_3_8,
    GeneralizedPauli_5_8,
    GeneralizedPauli_4_8,
    GeneralizedPauli_4_10,
    GeneralizedPauli_5_10,
    GeneralizedPauli_6_10,
    GeneralizedPauli_3_10,
    GeneralizedPauli_7_10
};

class v2RDMSolver: public Wavefunction{
  public:
    v2RDMSolver(SharedWavefunction reference_wavefunction,Options & options);
    ~v2RDMSolver();
    void common_init();

    double compute_energy();

    /// set constraints
    void set_constraints();

    /// number of primal variables (dimension of x)
    void determine_n_primal();

    /// number of constraints (dimension of y)
    void determine_n_dual();

    /// set offsets in x for each spin/symmetry block of rdms
    void set_primal_offsets();

    /// Au function for interfacing with bpsdp solver
    void bpsdp_Au(SharedVector A, SharedVector u);

    /// ATu function for interfacing with bpsdp solver
    void bpsdp_ATu(SharedVector A, SharedVector u);

    /// return dimension of the primal solution vector
    int n_primal(){return n_primal_;}

    /// return spin-free one-particle density matrix. full space. sparse
    std::vector<opdm> get_opdm_sparse(std::string type);

    /// return spin-free two-particle density matrix. full space. sparse.
    std::vector<tpdm> get_tpdm_sparse(std::string type);

    /// return spin-free one-particle density matrix as shared matrix for python API
    std::shared_ptr<Matrix> get_opdm();

    /// return spin-free two-particle density matrix as shared matrix for python API
    std::shared_ptr<Matrix> get_tpdm();

    /// return subset of orbitals for python API
    /**!
     * Similar to wavefunction.Ca_subset(); however, this version knows about all of the CI
     * subspaces in the SO basis. We stick to the definitions used by DETCI for now.
     * @param  orbital_name fzc, drc, docc, act, ras1, ras2, ras3, ras4, pop, vir, fzv, drv, or all
     * @return C            Returns the appropriate orbitals in the SO basis.
     */
    std::shared_ptr<Matrix> get_orbitals(const std::string &orbital_name);

    /**!
     * Similar to wavefunction.Ca_subset(); however, this version knows about all of the CI
     * subspaces in the SO basis. We stick to the MCSCF definitions for now.
     * @param  orbital_name FZC, DRC, DOCC, ACT, RAS1, RAS2, RAS3, RAS4, POP, VIR, FZV, DRV, or ALL
     * @param  orbitals     SharedMatrix to set
     */
    void set_orbitals(const std::string &orbital_name, SharedMatrix orbitals);

    /// Find out which orbitals belong where
    void orbital_locations(const std::string &orbital_name, int *start, int *end);

  protected:

    virtual bool same_a_b_orbs() const { return same_a_b_orbs_; }
    virtual bool same_a_b_dens() const { return same_a_b_dens_; }

    /// the sdp solver
    std::shared_ptr<BPSDPSolver> sdp_;

    double nalpha_;
    double nbeta_;

    /// constrain Q2 to be positive semidefinite?
    bool constrain_q2_;

    /// constrain G2 to be positive semidefinite?
    bool constrain_g2_;

    /// constraint T1 = D3 + Q3 to be positive semidefinite?
    bool constrain_t1_;

    /// constraint T2 = E3 + F3 to be positive semidefinite?
    bool constrain_t2_;

    /// keep d3 positive semidefinite and constrain D3->D2 mapping?
    bool constrain_d3_;

    /// keep d4 positive semidefinite and constrain D4->D3 mapping?
    bool constrain_d4_;

    /// constrain spin?
    bool constrain_spin_;

    /// constrain sz?
    bool constrain_sz_;

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
    long int n_dual_;

    /// total number of variables (dimension of primal solution vector)
    long int n_primal_;

    /// number of auxilliary basis functions
    long int nQ_;

    /// three-index integral buffer
    double * Qmo_;

    /// grab one-electron integrals (T+V) in MO basis
    SharedMatrix GetOEI();

    /// offsets
    int * d1aoff;
    int * d1boff;
    int * q1aoff;
    int * q1boff;
    int * d2aboff;
    int * d2aaoff;
    int * d2bboff;
    int * d200off;
    int * q2aboff;
    int * q2aaoff;
    int * q2bboff;
    int * g2aboff;
    int * g2baoff;
    int * g2aaoff;
    int * t1aaaoff;
    int * t1bbboff;
    int * t1aaboff;
    int * t1bbaoff;
    int * t2aaaoff;
    int * t2bbboff;
    int * t2aaboff;
    int * t2bbaoff;
    int * d3aaaoff;
    int * d3bbboff;
    int * d3aaboff;
    int * d3bbaoff;
    int * d4aaaaoff;
    int * d4aaaboff;
    int * d4aabboff;
    int * d4bbbaoff;
    int * d4bbbboff;

    /// standard vector of dimensions of each block of primal solution vector
    std::vector<int> dimensions_;

    /// standard vector of rank of each block of primal solution vector
    std::vector<int> rank_;

    int offset;

    // mapping arrays with abelian symmetry
    void build_mapping_arrays();
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

    /// triplets for each irrep:
    std::vector < std::vector < std::tuple<int,int,int> > > triplets;
    int * trip_aaa;
    int * trip_aab;
    int * trip_aba;
    int  ***  bas_aaa_sym;
    int  ***  bas_aab_sym;
    int  ***  bas_aba_sym;
    int **** ibas_aaa_sym;
    int **** ibas_aab_sym;
    int **** ibas_aba_sym;

    /// quartets for each irrep:
    std::vector < std::vector < std::tuple<int,int,int,int> > > quartets;
    int * quartet_aaaa;
    int * quartet_aaab;
    int * quartet_aabb;
    int  ***  bas_aaaa_sym;
    int  ***  bas_aaab_sym;
    int  ***  bas_aabb_sym;
    int ***** ibas_aaaa_sym;
    int ***** ibas_aaab_sym;
    int ***** ibas_aabb_sym;

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
    void T1_constraints_guess(SharedVector u);
    void T2_constraints_guess(SharedVector u);
    void Q2_constraints_guess(SharedVector u);
    void G2_constraints_guess(SharedVector u);

    void Spin_constraints_Au(SharedVector A,SharedVector u);
    void D2_constraints_Au(SharedVector A,SharedVector u);
    void Q2_constraints_Au(SharedVector A,SharedVector u);
    void G2_constraints_Au(SharedVector A,SharedVector u);
    void T1_constraints_Au(SharedVector A,SharedVector u);
    void T2_constraints_Au(SharedVector A,SharedVector u);
    void T2_constraints_Au_slow(SharedVector A,SharedVector u);
    void T2_tilde_constraints_Au(SharedVector A,SharedVector u);
    void D3_constraints_Au(SharedVector A,SharedVector u);
    void D4_constraints_Au(SharedVector A,SharedVector u);

    void Spin_constraints_ATu(SharedVector A,SharedVector u);
    void D2_constraints_ATu(SharedVector A,SharedVector u);
    void Q2_constraints_ATu(SharedVector A,SharedVector u);
    void G2_constraints_ATu(SharedVector A,SharedVector u);
    void T1_constraints_ATu(SharedVector A,SharedVector u);
    void T2_constraints_ATu(SharedVector A,SharedVector u);
    void T2_constraints_ATu_slow(SharedVector A,SharedVector u);
    void T2_tilde_constraints_ATu(SharedVector A,SharedVector u);
    void D3_constraints_ATu(SharedVector A,SharedVector u);
    void D4_constraints_ATu(SharedVector A,SharedVector u);


    // do enforce gpcs?
    bool constrain_gpc_;

    /// do enforce gpcs on 1rdm?
    bool constrain_gpc_1rdm_;

    /// do enforce gpcs on 2rdm?
    bool constrain_gpc_2rdm_;

    /// pick a set of gpcs based on na, nb, and amo
    void add_gpc_constraints(int na, int nb);

    /// generalized Pauli constraints
    void Generalized_Pauli_constraints_Au(SharedVector A,SharedVector u, int state);
    void Generalized_Pauli_3_8_constraints_Au(SharedVector A,SharedVector u, int state);
    void Generalized_Pauli_4_8_constraints_Au(SharedVector A,SharedVector u, int state);
    void Generalized_Pauli_3_6_constraints_Au(SharedVector A,SharedVector u, int state);
    void Generalized_Pauli_3_10_constraints_Au(SharedVector A,SharedVector u, int state);
    void Generalized_Pauli_4_10_constraints_Au(SharedVector A,SharedVector u, int state);
    void Generalized_Pauli_5_10_constraints_Au(SharedVector A,SharedVector u, int state);

    void Generalized_Pauli_constraints_ATu(SharedVector A,SharedVector u, int state);
    void Generalized_Pauli_3_8_constraints_ATu(SharedVector A,SharedVector u, int state);
    void Generalized_Pauli_4_8_constraints_ATu(SharedVector A,SharedVector u, int state);
    void Generalized_Pauli_3_6_constraints_ATu(SharedVector A,SharedVector u, int state);
    void Generalized_Pauli_3_10_constraints_ATu(SharedVector A,SharedVector u, int state);
    void Generalized_Pauli_4_10_constraints_ATu(SharedVector A,SharedVector u, int state);
    void Generalized_Pauli_5_10_constraints_ATu(SharedVector A,SharedVector u, int state);

    std::vector<std::shared_ptr<Matrix> > NatOrbs_;
    void SortedNaturalOrbitals(int state);

    double Generalized_Pauli_Au_term(double ** orbs,double * u,int * offa, int * offb,int index);
    void Generalized_Pauli_ATu_term(double val, double ** orbs,double * A,int * offa, int * offb,int index);

    double GP_N_8_Au(int & off, double * u, int * offa, int * offb, double ** orbs,
        double * eigvals, int d1, int d2, int d3, int d4, int d5, int d6, int d7,
        int d8);
    void GP_N_8_ATu(double dum,int & off, double * A, int * offa, int * offb,
        double ** orbs, int d1, int d2, int d3, int d4, int d5,
        int d6, int d7,int d8);

    double GP_N_10_Au(int & off, double * u, int * offa, int * offb, double ** orbs,
        double * eigvals, int d1, int d2, int d3, int d4, int d5, int d6, int d7,
        int d8, int d9, int d10);
    void GP_N_10_ATu(double dum,int & off, double * A, int * offa, int * offb,
        double ** orbs, int d1, int d2, int d3, int d4, int d5,
        int d6, int d7,int d8, int d9, int d10);

    /// offsets to gpc auxiliary variables in the primal vector
    std::vector<int * > gpcoff;

    /// the number of generalized pauli constraints
    int n_gpc_;

    /// the number of states to which gpcs are applied
    int n_gpc_states_; 

    /// which constraint will be applied
    std::vector<GeneralizedPauliConstraint> gpc_;

    /// do print gpc errors?
    bool print_gpc_error_;

// end of gpc

    /// SCF energy
    double escf_;

    /// nuclear repulsion energy
    double enuc_;

    //vectors
    SharedVector c;      // 1ei and 2ei of bpsdp
    SharedVector b;      // constraint vector
    SharedVector x;      // primal solution

    /// extended koopman's theorem
    void ExtendedKoopmans();
    void EKTEigensolver(std::shared_ptr<Matrix> V, std::shared_ptr<Matrix> D, std::shared_ptr<Vector> epsilon, bool use_dggev,std::string spin);

    /// compute natural orbitals and transform OPDM and TPDM to natural orbital basis
    void ComputeNaturalOrbitals();

    /// compute and print natural orbital occupation numbers
    void PrintNaturalOrbitalOccupations();

    /// push OPDM onto the wavefunction
    void FinalizeOPDM();

    // multiplicity
    int multiplicity_;

    /// full space of integrals for MO gradient / Hessian, blocked by symmetry
    double * tei_full_sym_;
    double * oei_full_sym_;
    // gidofalvi -- modified the type of tei_full_dim_ so that it is correct for large bases
    long int tei_full_dim_;
    int oei_full_dim_;
    std::shared_ptr<Matrix> T_;
    std::shared_ptr<Matrix> V_;

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
    void RepackIntegralsDF();

    /// compute frozen core energy and adjust oeis
    void FrozenCoreEnergy();

    /// function to rotate orbitals
    void RotateOrbitals();

    /// update ao/mo transformation matrix after orbital optimization
    // TODO: this function belongs in the new orbital optimizer class
    void UpdateTransformationMatrix();

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

    /// wall time for orbital optimization
    double orbopt_time_;

    /// total number of orbital optimization
    long int orbopt_iter_total_;

    /// write full 2RDM to disk
    void WriteTPDM();

    /// write full spin-free 2RDM to disk
    void WriteTPDMSpinFree();

    /// write full 1RDM to disk
    void WriteOPDM();

    /// write full 2RDM to disk in IWL format
    void WriteTPDM_IWL();

    /// write active-active-active-active 2RDM to disk
    void WriteActiveTPDM();

    /// read 2RDM from disk
    void ReadTPDM();

    /// write active 3RDM to disk
    void WriteActive3PDM();

    /// write molden file
    void WriteMoldenFile();

    /// read 3RDM from disk
    void Read3PDM();

    /// orbital lagrangian
    double * X_;

    void OrbitalLagrangian();

    /// memory available beyond what is allocated for v2RDM-CASSCF
    long int available_memory_;

    /// transform a four-index quantity from one basis to another
    void TransformFourIndex(double * inout, double * tmp, SharedMatrix trans);

    /// mo-mo transformation matrix
    SharedMatrix newMO_;

    /// FCIDUMP: dump integrals and RDMs to disk
    void FCIDUMP();

    /// break down energy into components
    void EnergyByComponent(double &kinetic, double &potential, double &two_electron_energy);

};

}
#endif

