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

#ifndef ORBITAL_OPTIMIZER_H
#define ORBITAL_OPTIMIZER_H

#include <string>

#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>

using namespace psi;

namespace hilbert{

class OrbitalOptimizer {

  public:

    OrbitalOptimizer(std::shared_ptr<Wavefunction> reference_wavefunction, Options & options);
    ~OrbitalOptimizer();

    /// take an orbital optimization step
    void optimize_orbitals(double * d2, double * d1, double * tei, double * oei, double * transformation_matrix,\
                           double& dE_orbopt, double& gnorm_orbopt, bool& converged_orbopt, int& iter_orbopt);

    /// returns the orbital lagrangian
    void get_lagrangian(std::shared_ptr<Matrix> Lagrangian);

    /// returns the orbital hessian
    void get_hessian(std::shared_ptr<Matrix> Hessian);

    /// has the optimization converged?
    bool is_converged(){ return is_energy_converged_ && is_gradient_converged_; }

  protected:

    /// is the energy converged to e_convergence_?
    bool is_energy_converged_;

    /// is the gradient converged to g_convergence_?
    bool is_gradient_converged_;

    /// number of irreps
    int nirrep_;

    /// number of auxiliary basis functions
    long int nQ_;

    /// number of molecular orbitals
    int nmo_;

    /// number of frozen core orbitals
    int nfrzc_;

    /// number of restricted core orbitals
    int nrstc_;

    /// number of active orbitals
    int amo_;

    /// number of restricted virtual orbitals
    int nrstv_;

    /// number of frozen virtual orbitals
    int nfrzv_;

    /// active orbitals per irrep
    int * amopi_;

    /// restricted core orbitals per irrep
    int * rstcpi_;

    /// restricted virtual orbitals per irrep
    int * rstvpi_;

    /// molecular orbitals per irrep
    Dimension nmopi_;

    /// doubly occuiped  orbitals per irrep
    Dimension doccpi_;

    /// frozen core orbitals per irrep
    Dimension frzcpi_;

    /// frozen virtual orbitals per irrep
    Dimension frzvpi_;

    /// orbital symmetries, in energy order
    int * orbital_symmetries_;

    /// energy convergence
    double e_convergence_;

    /// orbital gradient convergence
    double g_convergence_;

    /// maximum number of iterations
    int maxiter_;

    /// do print iteration information?
    bool write_;

    /// do evaluate the exact diagonal hessian?
    bool exact_diagonal_hessian_;

    /// do rotate active-active orbital pairs?
    bool active_active_rotations_;

    /// algorithm
    std::string algorithm_;

    /// variables GG added
 
    int * OindMap_e2i_;      // orbital index map energy--> irrep
    int * OindMap_e2c_;      // orbital index map energy --> class
    int * OindMap_c2e_;      // orbital index map class --> energy
    int * OindMap_c2i_;      // orbital index map class --> irrep
    int * OindMap_c2df_;     // orbital index map class --> df

    int * SymProd_;          // direct product table

    int * Active_GindMap_e_; // Geminal index (active orbitals indeces in energy order)
    int * Active_GindMap_c_; // Geminal index (active orbitals indeces in class order)
    int * Full_GindMap_c_;   // Geminal index (full orbitals indeces in class order)
    int * Grad_IndMap_c_;    // matrix to map orbital indeces (class order) on gradient orbital pair index

    int * act_gempi_;        // number of active geminals per irrep
    int * doc_gempi_;        // number of doubly occupied geminals per irrep
    int * ext_gempi_;        // number of external geminals per irrep
    int * full_gempi_;       // total number of geminals per irrep
    int * full_oe_offset_;   // offsets for indexing of symetry blocks of a 1-e matrix (full storage)
    int * lt_oe_offset_;     // offsets for indexing of symetry blocks of a 1-e matrix (LT storage)
    int * d2_irrep_offset_;  // offset array for indexing symmetry blocks of d2
    int * Q_offset_;         // offset array for indexing symmetry blocks of Q_
    int * U_offset_;         // offset array for indexing symmetry blocks of U_

    double * Fi_;               // inactive Fock matrix
    double * Fa_;               // active Fock matrix
    double * Q_;                // contraction of 2-e integrals with 2-e density 
    double * Z_;                // contraction of Fi with 1-e density
    double * orbital_gradient_; // orbital gradient
    double gradient_norm_;      // gradient norm
 
    double * kappa_current_;
    double * kappa_save_;
    double * kappa_ref_;
    double * kappa_new_;

    double * U_;             // transformation vectors
    double * Tei_scr1_;      // temporary array used in integral transform
    double * Tei_scr2_;      // temporary array used in integral transformation

    int ** first_index_;    
    int ** last_index_;

    int *** Fa_scr_offset_;

    double E1_c_;
    double E1_a_;
    double E2_cc_;
    double E2_ca_;
    double E2_aa_;
    double E_elec_;

    int Qstride_;
    int max_thread_ = 1;
    int nnz_Q_;
    int Nrot_;
 
    int U_dim_;
    int Tei_scr_dim_;

    double * diagonal_orbital_hessian_;

    // functions for optimization

    void Restore_last_step(double * tei, double * oei, double * transformation_matrix, bool converged_e, bool converged_g, bool eval_G);

    void Initiate_Optimizer(double * d2, double * d1, double * oei, double * transformation_matrix);

    void Finalize_Optimizer(double * d2, double * d1, double * oei, double * transformation_matrix);

    double Compute_EGH(double * d2, double * d1, double * tei, double * oei);

    bool Update_r(double dE_quad, double dE, double& r);

    void Update_kappa(double r_current);

    void Precondition_Step();

    double Compute_Quadratic_dE();

    void Alloc_Optimizer();
  
    void Dealloc_Optimizer();

    // functions for Hessian

    void Alloc_Hessian();

    void Dealloc_Hessian();

    void ComputeH_diag(double * d2, double * d1, double * tei );

    void ComputeH_diag_ad(double * d2, double * d1, double * tei );

    void tei_terms_H_ad(double * d2, double * d1, double * tei, int h);

    void ComputeH_diag_ed(double * tei); 

    void ComputeH_diag_aa(double * d2, double * d1, double * tei); 

    void ComputeH_diag_ea(double * d2, double * d1, double * tei);  

    void tei_terms_H_ea(double * d2, double * tei, int h);

    void PrintH_diag();

    void CleanH_diag();

    // functions to determine U = e^(-kappa)

    void ExponentiateKappa(double * kappa);

    void ExponentiateKappa_UnpackSymBlock(double * kappa, int h);

    void ExponentiateKappa_ComputeUSymBlock(int h);

    // Funfctions for integral transformation

    void TransformIntegrals(double * tei, double * oei, double * MOcoeff);

    void TransformOei(double * oei);

    void TransformMOcoeff(double * MOcoeff);

    void Print_IntQ(double * tei);

    void TransformOei_UnpackSymBlock(double * oei, int h );

    void TransformOei_PackSymBlock(double * oei, int h );

    void TransformTei(double * tei);

    void TransformTei_Q(double * teiQ, double * tei_scr1, double * tei_scr2);

    void TransformTei_UnpackQ(double * teiQ, double * tei_scr1, int h_p, int h_q);

    void TransformTei_PackQ(double * teiQ, double * tei_scr1, int h_p, int h_q);

    // functions for calculating gradient 
 
    void ComputeG(double * d2, double * d1, double * tei, double * oei);
    
    void Compute_Fi(double * tei, double * oei);
   
    void Compute_Fi_J(double * tei, double * oei);

    void Compute_Fi_K(double * tei);

    void Compute_Fa(double * tei, double * d1);

    void Compute_Fa_J(double * tei, double * d1);

    void Compute_Fa_K(double * tei, double * d1);

    void Compute_Fa_K_copy_3index(double * tei, double * scr_3index_1, \
                                  int h_p, int p_class);

    void Compute_Fa_K_d_block(double * scr_4index_1, double * scr_3index_1, \
                              double * d1, int h_p, int p_class, int max_dim);

    void Compute_Fa_K_od_block(double * scr_4index_1, double * scr_3index_1, \
                               double * scr_3index_2, double * d1, int h_p,  \
                               int p_class, int q_class, int max_dim);

    void Alloc_Fa_scr_offset();

    void Dealloc_Fa_scr_offset();

    void Compute_Fa_scr_offset();

    void Compute_Q(double * tei, double * d2);

    void Initialize_scratch_G();

    void PrintQ(double * Q);

    void Compute_Q_update(double * tei, double * scr_3index_1, \
                          double * scr_3index_2, double * scr_4index_1, int h_tu);

    void Compute_Q_copy_d2(double * d2, double * scr_4index_2, int h_tu);

    void Compute_Q_copy_tei(double * tei, double * scr_3index_1, int h_tu);

    void Compute_Z(double * d1);

    void Compute_G();

    void Compute_G_ed();

    void Compute_G_ea();

    void Compute_G_aa();

    void Compute_G_ad();

    void Determine_GradIndMap();

    void Printfock(double * fock);

    // functions for calculating energy

    double ComputeE_2index(double * d1, double * oei);

    double ComputeE(double * d2, double * d1, double * tei, double * oei);

    double ComputeE2_aa(double * d2, double * tei);

    double ComputeE2_ca(double * d1, double * tei);

    double ComputeE2_cc(double * tei);

    double ComputeE1_c(double * oei);

    double ComputeE1_a(double * d1, double * oei);

    // functions for allocating and determining variables and mapping arrays

    void Determine_OindMap();

    void Alloc_OindMap();

    void Dealloc_OindMap();

    void Alloc_Gradient();

    void Dealloc_Gradient();

    void Alloc_U();

    void Dealloc_U();

    // various helper functions

    int GemSym(int h_i, int h_j);

    int Full_ij_index(int row, int col, int nrow);

    int Lt_ij_index(int row, int col);

    int Full_ij_oe_index( int p_c, int q_c, int h);

    int Lt_ij_oe_index( int p_c, int q_c, int h);

    int df_ij_index( int p_c, int q_c);

    int Q_index( int t_c, int p_c, int h_p );
 
    int Max_value( int * array, int dim);

    // functions for sorting & scaling d2, d1, and e1

    void Sort(double * d2, double * d1, double * oei, double * transformation_matrix, int direction);

};

}

#endif
