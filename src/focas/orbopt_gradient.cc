#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsi4util/PsiOutStream.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/molecule.h>
#include <psi4/libmints/basisset.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/psifiles.h>
#include <psi4/libqt/qt.h>

#include <misc/blas.h>
#include <misc/omp.h>

#include "orbital_optimizer.h"

using namespace psi;
using namespace fnocc;

namespace hilbert{


void OrbitalOptimizer::ComputeG(double * d2, double * d1, double * tei, double * oei){

    Initialize_scratch_G();

    Compute_Fi(tei, oei);

    Compute_Fa(tei, d1);

    Compute_Q(tei, d2);

    Compute_Z(d1);

    Compute_G();

//    if ( RAS_aa_ ) RAS_ZeroRedundant_aa(orbital_gradient_);

    gradient_norm_ = sqrt( C_DDOT(Nrot_,&orbital_gradient_[0],1,&orbital_gradient_[0],1) );

}

void OrbitalOptimizer::Compute_G(){

    Compute_G_ea();

    Compute_G_ed();

    if ( active_active_rotations_ ) Compute_G_aa();

    Compute_G_ad();

}

void OrbitalOptimizer:: Compute_G_ed(){

    for ( int h = 0; h < nirrep_; h++ ){

        for ( int a_c = first_index_[h][2]; a_c < last_index_[h][2]; a_c++ ){

            for ( int i_c = first_index_[h][0]; i_c < last_index_[h][0]; i_c++ ){

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(a_c,i_c)];

                orbital_gradient_[grad_ind] = 4.0e0 * ( Fi_[Lt_ij_oe_index(a_c,i_c,h) ] + \
                                                        Fa_[Lt_ij_oe_index(a_c,i_c,h) ] );

            } // end q_c loop 

        } // end p_c loop

    } // end h loop

}

void OrbitalOptimizer:: Compute_G_ea(){

    for ( int h = 0; h < nirrep_; h++ ){

        for ( int a_c = first_index_[h][2]; a_c < last_index_[h][2]; a_c++ ){

            for ( int t_c = first_index_[h][1]; t_c < last_index_[h][1]; t_c++ ){

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(a_c,t_c)];

                orbital_gradient_[grad_ind] = 2.0e0 * ( Q_[Q_index(t_c,a_c,h)] + \
                                                        Z_[Q_index(t_c,a_c,h)] );

            } // end q_c loop 

        } // end p_c loop

    } // end h loop

}

void OrbitalOptimizer:: Compute_G_aa(){

    for ( int h = 0; h < nirrep_; h++ ){

        for ( int t_c = first_index_[h][1]; t_c < last_index_[h][1]; t_c++ ){

            for ( int u_c = first_index_[h][1]; u_c < t_c; u_c++ ){

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(t_c,u_c)];

                orbital_gradient_[grad_ind] = 2.0e0 * ( Q_[Q_index(u_c,t_c,h)] + \
                                                        Z_[Q_index(u_c,t_c,h)] - \
                                                        Q_[Q_index(t_c,u_c,h)] - \
                                                        Z_[Q_index(t_c,u_c,h)] );

            } // end q_c loop 

        } // end p_c loop

    } // end h loop

}

void OrbitalOptimizer:: Compute_G_ad(){

    for ( int h = 0; h < nirrep_; h++ ){

        for ( int t_c = first_index_[h][1]; t_c < last_index_[h][1]; t_c++ ){

            for ( int i_c = first_index_[h][0]; i_c < last_index_[h][0]; i_c++ ){

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(t_c,i_c)];

                orbital_gradient_[grad_ind] = 4.0e0 * ( Fi_[Lt_ij_oe_index(t_c,i_c,h) ] +   \
                                                        Fa_[Lt_ij_oe_index(t_c,i_c,h) ] ) - \
                                              2.0e0 * ( Q_[Q_index(t_c,i_c,h)] + \
                                                        Z_[Q_index(t_c,i_c,h)] );

            } // end q_c loop 

        } // end p_c loop

    } // end h loop

}


void OrbitalOptimizer::Determine_GradIndMap(){

    int Nrot = 0;

    for ( int p_class = 0; p_class < 3; p_class++ ){

        int q_class_max = p_class;
        if ( active_active_rotations_ && p_class == 1 ) q_class_max = p_class + 1;

        for ( int q_class = 0; q_class < q_class_max; q_class++ ){

            for ( int h = 0; h < nirrep_; h++ ){

                for ( int p_c = first_index_[h][p_class]; p_c < last_index_[h][p_class]; p_c++ ){

                    int q_c_max = last_index_[h][q_class];
                    if ( p_class == q_class ) q_c_max = p_c;

                    for ( int q_c = first_index_[h][q_class]; q_c < q_c_max; q_c++ ){

                        Grad_IndMap_c_[Lt_ij_index(p_c,q_c)] = Nrot;

                        Nrot++;

                    } // end q_c loop 

                } // end p_c loop

            } // end h loop

        } // emd q_class loop       

    } // end p_class loop

    if ( Nrot != Nrot_ ) printf("!!! Error ... Nrot mismatch !!! Nrot=%i Nrot_=%i\n",Nrot,Nrot_);

}

void OrbitalOptimizer::Compute_Z(double * d1){

    for ( int h_t = 0; h_t < nirrep_; h_t++){

        for ( int t_c = first_index_[h_t][1]; t_c < last_index_[h_t][1]; t_c++){

            for ( int p_class =0; p_class < 3; p_class++){

                for ( int p_c = first_index_[h_t][p_class]; p_c < last_index_[h_t][p_class]; p_c++ ){

                    double Zval = 0.0e0;

                    for ( int u_c = first_index_[h_t][1]; u_c < last_index_[h_t][1]; u_c++ ){

                        int tu_den = Active_GindMap_c_[Full_ij_index(t_c,u_c,nmo_)];
                        int pu_fi  = Lt_ij_oe_index(p_c,u_c,h_t);

                        Zval      += d1[tu_den] * Fi_[pu_fi];

                    } // end u_c loop 

                    Z_[Q_index(t_c,p_c,h_t)] += Zval;

                } // end p_c loop

            } // end p_class loop

        } // end t_c loop

    } // end h_p loop

}

void OrbitalOptimizer::Compute_Q(double * tei, double * d2){

    int max_gempi = Max_value(act_gempi_,nirrep_);
    int max_act   = Max_value(amopi_,nirrep_);

    double * scr_4index_1 = (double *)malloc(max_gempi*max_gempi*sizeof(double));
    memset((void*)scr_4index_1,'\0',max_gempi*max_gempi*sizeof(double));

    double * scr_3index_2 = (double *)malloc(max_gempi*nQ_*sizeof(double));
    memset((void*)scr_3index_2,'\0',max_gempi*nQ_*sizeof(double));

    double * scr_3index_1 = (double *)malloc(max_gempi*nQ_*sizeof(double));
    memset((void*)scr_3index_1,'\0',max_gempi*nQ_*sizeof(double));

    for ( int h_tu = 0; h_tu < nirrep_; h_tu++){

        int ngem_tu = act_gempi_[h_tu];

        if ( ngem_tu == 0 ) continue;

        // unpack density and scale off-diagonal elements

        Compute_Q_copy_d2(d2, scr_4index_1, h_tu);

        // copy 3-index integrals

        Compute_Q_copy_tei(tei, scr_3index_1, h_tu);

        // contract integrals with density

        F_DGEMM('n','n',nQ_,ngem_tu,ngem_tu,1.0e0,&scr_3index_1[0],nQ_,\
                                                  &scr_4index_1[0],ngem_tu,\
                                            0.0e0,&scr_3index_2[0],nQ_);

        // Update Q matrix

        Compute_Q_update(tei, scr_3index_1, scr_3index_2, scr_4index_1, h_tu);

    }

    free(scr_4index_1);
    free(scr_3index_1);
    free(scr_3index_2);

}

void OrbitalOptimizer::PrintQ(double * Q){

    for ( int h_p = 0; h_p < nirrep_; h_p++){

        for ( int p_class = 0; p_class < 3; p_class++){

            for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                for ( int t_c = first_index_[h_p][1]; t_c < last_index_[h_p][1]; t_c++){

                    printf("%i %i   %20.12e\n",t_c,p_c,Q[Q_index(t_c,p_c,h_p)]);

                }

            }

        }

    }

}

void OrbitalOptimizer::Compute_Q_update(double * tei, double * scr_3index_1, double * scr_3index_2, double * scr_4index_1, int h_tu){

    int ngem_tu = act_gempi_[h_tu];

    if ( ngem_tu == 0 ) return;

    for ( int h_p = 0; h_p < nirrep_; h_p++){

        int h_u   = GemSym(h_p,h_tu);

        int amo_u = amopi_[h_u];

        if ( amo_u == 0 || nmopi_[h_p] == 0 ) continue;

        for (int p_class = 0; p_class < 3; p_class++){

            for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                for ( int u_c = first_index_[h_u][1]; u_c < last_index_[h_u][1]; u_c++){

                    int pu_df   = df_ij_index(p_c,u_c);

                    int scr_ind = ( u_c - first_index_[h_u][1] ) * nQ_;

                    C_DCOPY(nQ_,&tei[pu_df],Qstride_,&scr_3index_1[scr_ind],1);

                } // end u_c loop

                F_DGEMM('t','n',amo_u,ngem_tu,nQ_,1.0e0,&scr_3index_1[0],nQ_,\
                                                        &scr_3index_2[0],nQ_,\
                                                  0.0e0,&scr_4index_1[0],amo_u);

                for ( int t_c = first_index_[h_p][1]; t_c < last_index_[h_p][1]; t_c++){

                    double Qval = 0.0e0;

                    int u_scr = 0;

                    for ( int u_c= first_index_[h_u][1]; u_c < last_index_[h_u][1]; u_c++){

                        int tu_den = Active_GindMap_c_[Full_ij_index(t_c,u_c,nmo_)];

                        Qval += scr_4index_1[Full_ij_index(u_scr,tu_den,amo_u)];

                        u_scr++;

                    } // end u_c loop

                    Q_[Q_index(t_c,p_c,h_p)] += Qval;

                } // end t_c loop

            } // end p_c loop

        } // end p_class loop

    } // end h_p loop

}

void OrbitalOptimizer::Compute_Q_copy_tei(double * tei, double * scr_3index_1, int h_tu){

    int ngem_tu = act_gempi_[h_tu];

    if ( ngem_tu == 0 ) return;

    for ( int h_v = 0; h_v < nirrep_; h_v++){

        int h_w = GemSym(h_v,h_tu);

        if ( h_w > h_v ) continue;

        for ( int v_c = first_index_[h_v][1]; v_c < last_index_[h_v][1]; v_c++){

            int max_w_c = last_index_[h_w][1];
            if ( h_tu == 0 ) max_w_c = v_c + 1;

            for ( int w_c = first_index_[h_w][1]; w_c < max_w_c; w_c++){

                int vw_df  = df_ij_index(w_c,v_c);
                int vw_scr = Active_GindMap_c_[Full_ij_index(v_c,w_c,nmo_)] * nQ_;

                C_DCOPY(nQ_,&tei[vw_df],Qstride_,&scr_3index_1[vw_scr],1);

            } // end w_c loop

        } // end v_c loop

    } // end h_v loop

}

void OrbitalOptimizer::Compute_Q_copy_d2(double * d2, double * scr_4index_2, int h_tu){

    int ngem_tu = act_gempi_[h_tu];

    int tu_offset = d2_irrep_offset_[h_tu];

    if ( ngem_tu == 0 ) return;

    for ( int h_t = 0; h_t < nirrep_; h_t++){

        int h_u = GemSym(h_t,h_tu);

        if ( h_u > h_t ) continue;

        for ( int t_c = first_index_[h_t][1]; t_c < last_index_[h_t][1]; t_c++){

            int max_u_c = last_index_[h_u][1];
            if ( h_tu == 0 ) max_u_c = t_c + 1;

            for ( int u_c = first_index_[h_u][1]; u_c < max_u_c; u_c++){

                int tu_den = Active_GindMap_c_[Full_ij_index(t_c,u_c,nmo_)];

                for ( int h_v = 0; h_v < nirrep_; h_v++){

                    int h_w = GemSym(h_v,h_tu);

                    if ( h_w > h_v ) continue;

                    for ( int v_c = first_index_[h_v][1]; v_c < last_index_[h_v][1]; v_c++){

                        int max_w_c = last_index_[h_w][1];
                        if ( h_tu == 0 ) max_w_c = v_c;

                        for ( int w_c = first_index_[h_w][1]; w_c < max_w_c; w_c++){

                            int vw_den   = Active_GindMap_c_[Full_ij_index(v_c,w_c,nmo_)];
                            int vwtu_den = Lt_ij_index(tu_den,vw_den) + tu_offset;

                            scr_4index_2[Full_ij_index(vw_den,tu_den,ngem_tu)] = 2.0e0*d2[vwtu_den];

                        } // end w_c loop

                        if ( h_tu != 0 ) continue;

                        int vv_den = Active_GindMap_c_[Full_ij_index(v_c,v_c,nmo_)];
                        int vvtu_den = Lt_ij_index(vv_den,tu_den) + tu_offset;

                        scr_4index_2[Full_ij_index(vv_den,tu_den,ngem_tu)] = d2[vvtu_den];

                    } // end v_c loop

                } // end h_v loop     

            } // end u_c loop

        } // end t_c loop

    } // end h_t loop

}


void OrbitalOptimizer::Compute_Fa(double * tei, double * d1){

    Compute_Fa_J(tei, d1);

    Compute_Fa_K(tei, d1);

}

void OrbitalOptimizer::Compute_Fa_K(double * tei, double *d1){

    // function to calculate exchange contribution to Fa

    int max_act = Max_value(amopi_,nirrep_);
    int max_doc = Max_value(rstcpi_,nirrep_);
    int max_ext = Max_value(rstvpi_,nirrep_);

    int max_dim = max_doc; if ( max_act > max_dim) max_dim = max_act; if ( max_ext > max_dim) max_dim = max_ext;

    int max_dim_4ind = max_dim;

    max_dim *= amo_;

    double * scr_3index_1 = (double *)malloc(max_dim*nQ_*sizeof(double));
    memset((void*)scr_3index_1,'\0',max_dim*nQ_*sizeof(double));

    double * scr_3index_2 = (double *)malloc(max_dim*nQ_*sizeof(double));
    memset((void*)scr_3index_2,'\0',max_dim*nQ_*sizeof(double));

    double * scr_4index_1 = (double *)malloc(max_dim_4ind*max_dim_4ind*sizeof(double));
    memset((void*)scr_4index_1,'\0',max_dim_4ind*max_dim_4ind*sizeof(double));

    Alloc_Fa_scr_offset();

    Compute_Fa_scr_offset();

    for ( int h_p = 0; h_p < nirrep_; h_p++){

        for ( int p_class = 0; p_class < 3; p_class++){

            int nmo_p = last_index_[h_p][p_class] - first_index_[h_p][p_class];

            if ( nmo_p == 0 ) continue;

            // gather all (:,pt) for all p in this h_p and p_class & all active t

            Compute_Fa_K_copy_3index(tei, scr_3index_1, h_p, p_class);

            // Calculate diagonal elements of F_a

            Compute_Fa_K_d_block(scr_4index_1, scr_3index_1, d1, h_p, p_class, max_dim_4ind);

            // Calculate off-diagonal element p_class & q_class = p_class -1 

            int q_class = p_class - 1;

            if ( p_class == 1 || p_class == 2 ) Compute_Fa_K_od_block(scr_4index_1, scr_3index_1, scr_3index_2, d1, h_p, p_class, q_class, max_dim_4ind); 

            // save these 3index integrals

            if ( p_class < 2 ) C_DCOPY(nmo_p*amo_*nQ_,&scr_3index_1[0],1,&scr_3index_2[0],1);

if ( p_class == 2 ) {

                // at this point scr_3index_1 has p_class == 2 ; gather q_class = 0 into scr_index_2

                q_class = 0;

                Compute_Fa_K_copy_3index(tei, scr_3index_2, h_p, q_class);

                Compute_Fa_K_od_block(scr_4index_1, scr_3index_1, scr_3index_2, d1, h_p, p_class, q_class, max_dim_4ind);
            }

        } // end p_class loop

    } // end h_p loop

    Dealloc_Fa_scr_offset();

    free(scr_3index_1);
    free(scr_3index_2);
    free(scr_4index_1);
}

void OrbitalOptimizer::Compute_Fa_K_copy_3index(double * tei, double * scr_3index_1, int h_p, int p_class){

    int nmo_p = last_index_[h_p][p_class] - first_index_[h_p][p_class];

    if ( nmo_p == 0 ) return;

    for ( int h_t = 0; h_t < nirrep_; h_t++){

        int pt_scr_offset_1 = Fa_scr_offset_[h_p][p_class][h_t];

        for ( int t_c = first_index_[h_t][1]; t_c < last_index_[h_t][1]; t_c++){

            int pt_scr_offset_2 = pt_scr_offset_1 + ( t_c - first_index_[h_t][1] ) * nmo_p;

            for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                int pt_df  = df_ij_index(p_c,t_c);

                int pt_scr = ( pt_scr_offset_2 + p_c - first_index_[h_p][p_class] ) * nQ_;

                C_DCOPY(nQ_,&tei[pt_df],Qstride_,&scr_3index_1[pt_scr],1);

           } // end p_c loop

        } // end t_c loop

    } // end h_t loop

}

void OrbitalOptimizer::Compute_Fa_K_d_block(double * scr_4index_1, double * scr_3index_1, \
                                             double * d1, int h_p, int p_class, int max_dim){

    int nmo_p = last_index_[h_p][p_class] - first_index_[h_p][p_class];

    if ( nmo_p == 0 ) return;

    if ( p_class == 0 || p_class == 1 ) {

        memset((void*)scr_4index_1,'\0',max_dim*max_dim*sizeof(double));

        for ( int h_t = 0; h_t < nirrep_; h_t++){

            int nmo_t = last_index_[h_t][1] - first_index_[h_t][1];

            int pt_scr_offset = Fa_scr_offset_[h_p][p_class][h_t];
            int qu_scr_offset = pt_scr_offset;

            for ( int t_c = first_index_[h_t][1]; t_c < last_index_[h_t][1]; t_c++){

                int pt_scr_ind = ( pt_scr_offset + ( t_c - first_index_[h_t][1] ) * nmo_p ) * nQ_;

                for ( int u_c = first_index_[h_t][1]; u_c < last_index_[h_t][1]; u_c++){

                    int qu_scr_ind = ( qu_scr_offset + ( u_c - first_index_[h_t][1] ) * nmo_p ) * nQ_;
                    int tu_den     = Active_GindMap_c_[Full_ij_index(t_c,u_c,nmo_)];

                    double d1_val = 5.0e-1*d1[tu_den];

                    F_DGEMM('t','n',nmo_p,nmo_p,nQ_,d1_val,&scr_3index_1[pt_scr_ind],nQ_,\
                                                           &scr_3index_1[qu_scr_ind],nQ_,\
                                                           1.0e0,&scr_4index_1[0],nmo_p);

                } // end u_c loop 

            } // end t_c loop

        } // end h_t loop

        for ( int q_c = first_index_[h_p][p_class]; q_c < last_index_[h_p][p_class]; q_c++){

            int q_scr = q_c - first_index_[h_p][p_class];

            for ( int p_c = q_c; p_c < last_index_[h_p][p_class]; p_c++){

                int p_scr    = p_c - first_index_[h_p][p_class];
                int Fa_index = Lt_ij_oe_index(p_c,q_c,h_p);

                Fa_[Fa_index] -= scr_4index_1[Full_ij_index(p_scr,q_scr,nmo_p)];

            } // end q_c loop

        } // end q_c loop

    }

    if ( p_class == 2) {

        for ( int h_t = 0; h_t < nirrep_; h_t++){

            int scr_offset = Fa_scr_offset_[h_p][p_class][h_t];

            for ( int t_c = first_index_[h_t][1]; t_c < last_index_[h_t][1]; t_c++){

                int pt_scr_offset = scr_offset + ( t_c - first_index_[h_t][1] ) * nmo_p;

                for ( int u_c = first_index_[h_t][1]; u_c < last_index_[h_t][1]; u_c++){

                    int pu_scr_offset = scr_offset + ( u_c - first_index_[h_t][1] ) * nmo_p;
                    int tu_den        = Active_GindMap_c_[Full_ij_index(t_c,u_c,nmo_)];

                    double d1_val = 5.0e-1*d1[tu_den];

                    int p_red = 0;

                    for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                        int pt_scr_index = ( pt_scr_offset + p_red ) * nQ_;
                        int pu_scr_index = ( pu_scr_offset + p_red ) * nQ_;
                        int Fa_index     = Lt_ij_oe_index(p_c,p_c,h_p);

                        Fa_[Fa_index] -= d1_val * C_DDOT(nQ_,&scr_3index_1[pt_scr_index],1,&scr_3index_1[pu_scr_index],1);
                        p_red++;


                    } // end p_c loop

                } // end u_c loop

            } // end t_c loop

        } // end h_t loop  

    }

}

void OrbitalOptimizer::Compute_Fa_K_od_block(double * scr_4index_1, double * scr_3index_1, \
                                             double * scr_3index_2, double * d1, int h_p,  \
                                             int p_class, int q_class, int max_dim){

    memset((void*)scr_4index_1,'\0',max_dim*max_dim*sizeof(double));

    int nmo_p = last_index_[h_p][p_class] - first_index_[h_p][p_class];
    int nmo_q = last_index_[h_p][q_class] - first_index_[h_p][q_class];

    if ( nmo_q == 0 || nmo_q == 0 ) return;

    for ( int h_t = 0; h_t < nirrep_; h_t++){

        int nmo_t = last_index_[h_t][1] - first_index_[h_t][1];

        if ( nmo_t == 0 ) continue;

        int pt_scr_offset = Fa_scr_offset_[h_p][p_class][h_t];
        int qu_scr_offset = Fa_scr_offset_[h_p][q_class][h_t];

        for ( int u_c = first_index_[h_t][1]; u_c < last_index_[h_t][1]; u_c++){

            int qu_scr_ind = ( qu_scr_offset + ( u_c - first_index_[h_t][1] ) * nmo_q ) * nQ_;

            int qu_t = qu_scr_ind / nQ_;

            for ( int t_c = first_index_[h_t][1]; t_c < last_index_[h_t][1]; t_c++){

                 int tu_den        = Active_GindMap_c_[Full_ij_index(t_c,u_c,nmo_)];

                 double d1_val = 5.0e-1*d1[tu_den];

                 int pt_scr_ind = ( pt_scr_offset + ( t_c - first_index_[h_t][1] ) * nmo_p ) * nQ_;

                 int pt_t = pt_scr_ind / nQ_;

                 F_DGEMM('t','n',nmo_p,nmo_q,nQ_,d1_val,&scr_3index_1[pt_scr_ind],nQ_,\
                                                         &scr_3index_2[qu_scr_ind],nQ_,\
                                                          1.0e0,&scr_4index_1[0],nmo_p);

            } // end u_c loop

        } // end t_c loop    

    } // end h_t loop

    for ( int q_c = first_index_[h_p][q_class]; q_c < last_index_[h_p][q_class]; q_c++){

         int q_scr = q_c - first_index_[h_p][q_class];

         for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

             int p_scr    = p_c - first_index_[h_p][p_class];
             int Fa_index = Lt_ij_oe_index(p_c,q_c,h_p);

             Fa_[Fa_index] -= scr_4index_1[Full_ij_index(p_scr,q_scr,nmo_p)];

         } // end p_c loop

    } // end q_c loop

}

void OrbitalOptimizer::Alloc_Fa_scr_offset(){

    Fa_scr_offset_ = (int***)malloc(nirrep_*sizeof(int**));

    for (int h_p = 0; h_p < nirrep_; h_p++){

        Fa_scr_offset_[h_p] = (int**)malloc(3*sizeof(int*));

        for ( int p_class = 0; p_class < 3; p_class++){

            Fa_scr_offset_[h_p][p_class] = (int*)malloc(nirrep_*sizeof(int));

        }
    }

}

void OrbitalOptimizer::Dealloc_Fa_scr_offset(){

    for ( int h_p = 0; h_p < nirrep_; h_p++){

        for ( int p_class = 0; p_class < 3; p_class++){

            free(Fa_scr_offset_[h_p][p_class]);

        }

        free(Fa_scr_offset_[h_p]);

    }

    free(Fa_scr_offset_);

}

void OrbitalOptimizer::Compute_Fa_scr_offset(){

    for (int h_p = 0; h_p < nirrep_; h_p++){

        for ( int p_class = 0; p_class < 3; p_class++){

            int dim = 0;

            int nmo_p = rstcpi_[h_p];
            if ( p_class == 1 ) nmo_p = amopi_[h_p];
            if ( p_class == 2 ) nmo_p = rstvpi_[h_p];

            for ( int h_t = 0; h_t < nirrep_; h_t++){

                Fa_scr_offset_[h_p][p_class][h_t] = dim;

                dim += nmo_p * amopi_[h_t];

            }

        }
    }

}

void OrbitalOptimizer::Compute_Fa_J(double * tei, double *d1){

    // function to calculate Coulomb contribution to Fa 
    // set up for parallel execution (thread_id should be determined with openmp calls)

    double * scr_3index_1 = (double *)malloc(nQ_*sizeof(double));
    memset((void*)scr_3index_1,'\0',nQ_*sizeof(double));

    double * scr_3index_2 = (double *)malloc(max_thread_*nQ_*sizeof(double));
    memset((void*)scr_3index_2,'\0',max_thread_*nQ_*sizeof(double));

    for ( int h_t = 0; h_t < nirrep_; h_t++){

        for ( int t_c = first_index_[h_t][1]; t_c < last_index_[h_t][1]; t_c++){

            int thread_id = 0;

            for ( int u_c = first_index_[h_t][1]; u_c < t_c; u_c++){

                int tu_df  = df_ij_index(t_c,u_c);
                int tu_den = Active_GindMap_c_[Full_ij_index(t_c,u_c,nmo_)];

                C_DAXPY(nQ_,2.0e0*d1[tu_den],&tei[tu_df],Qstride_,&scr_3index_2[thread_id*nQ_],1);

            } // end u_c loop

            int tt_df  = df_ij_index(t_c,t_c);
            int tt_den = Active_GindMap_c_[Full_ij_index(t_c,t_c,nmo_)];

            C_DAXPY(nQ_,d1[tt_den],&tei[tt_df],Qstride_,&scr_3index_2[thread_id*nQ_],1);

        } // end t_c loop

    } // end h_t loop

    for ( int thread_reduce = 0; thread_reduce < max_thread_;thread_reduce++){

        C_DAXPY(nQ_,1.0e0,&scr_3index_2[thread_reduce*nQ_],1,&scr_3index_1[0],1);

    } // end thread_reduce loop

    for ( int p_class = 0; p_class < 3; p_class++){

        for ( int q_class = 0; q_class <= p_class; q_class++){

            for ( int h_p = 0; h_p < nirrep_; h_p++){

                for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                    int q_c_min = first_index_[h_p][q_class];
                    if ( q_class == 2 ) q_c_min = p_c;

                    int q_c_max = last_index_[h_p][q_class];
                    if ( q_class == p_class ) q_c_max = p_c + 1;

                    for ( int q_c = q_c_min; q_c < q_c_max; q_c++){

                        int pq_df    = df_ij_index(p_c,q_c);
                        int Fa_index = Lt_ij_oe_index(p_c,q_c,h_p);

                        Fa_[Fa_index] = C_DDOT(nQ_,&tei[pq_df],Qstride_,&scr_3index_1[0],1);

                    } // end q_c loop

                } // end p_c loop

            } // end h_p loop

        } // end q_class loop

    } // end p_class loop

    free(scr_3index_1);
    free(scr_3index_2);

}


void OrbitalOptimizer::Compute_Fi(double * tei, double * oei){

    Compute_Fi_J(tei, oei);

    Compute_Fi_K(tei);

}

void OrbitalOptimizer::Compute_Fi_K(double * tei){

    // function to calculate Exchange contribution to Fi
    // set up so that it can run in parallel (as long as thread_id is called with openmp function)
    // currently, only diagonal elements are calculate for e-e block (can be changed by the q_c_min variable)

    int max_doc = Max_value(rstcpi_,nirrep_);
    int max_act = Max_value(amopi_,nirrep_);
    int max_ext = Max_value(rstvpi_,nirrep_);

    int max_dim = max_doc; if ( max_act > max_dim) max_dim = max_act; if ( max_ext > max_dim) max_dim = max_ext;

    double * scr_3index_1 = (double *)malloc(max_dim*nQ_*sizeof(double));
    memset((void*)scr_3index_1,'\0',max_dim*nQ_*sizeof(double));

    double * scr_3index_2 = (double *)malloc(max_dim*nQ_*sizeof(double));
    memset((void*)scr_3index_2,'\0',max_dim*nQ_*sizeof(double));

    double * scr_4index_1 = (double *)malloc(max_dim*max_dim*sizeof(double));
    memset((void*)scr_3index_1,'\0',max_dim*max_dim*sizeof(double));

    double * scr_4index_2 = (double *)malloc(max_dim*max_dim*sizeof(double));
    memset((void*)scr_3index_2,'\0',max_dim*max_dim*sizeof(double));

    double * scr_diag_ext = (double *)malloc(max_dim*max_thread_*sizeof(double));
    memset((void*)scr_diag_ext,'\0',max_dim*max_thread_*sizeof(double));

    int * calc_p_diag = (int *)malloc(3*nirrep_*sizeof(int));

    int * calc_q_diag = (int *)malloc(3*nirrep_*sizeof(int));

    for ( int p=0;p<3;p++){

        for ( int h=0;h<nirrep_;h++){

            calc_p_diag[p*nirrep_ + h] = 0;
            calc_q_diag[p*nirrep_ + h] = 1;

        }
    }

    for ( int h = 0; h < nirrep_; h++) calc_q_diag[2*nirrep_+h] = 0;
    for ( int h = 0; h < nirrep_; h++) calc_p_diag[2*nirrep_+h] = 1;

    for ( int p_class = 0; p_class < 3; p_class++){

        for ( int q_class = 0; q_class < p_class; q_class++){

            for ( int h_p = 0; h_p < nirrep_; h_p++){

                int nmo_p = last_index_[h_p][p_class] - first_index_[h_p][p_class];

                int nmo_q = last_index_[h_p][q_class] - first_index_[h_p][q_class];

                memset((void*)scr_4index_1,'\0',max_dim*max_dim*sizeof(double));
                if ( calc_q_diag[q_class*nirrep_ + h_p] == 1 ) memset((void*)scr_4index_2,'\0',max_dim*max_dim*sizeof(double));
                if ( calc_p_diag[p_class*nirrep_ + h_p] == 1 ) memset((void*)scr_diag_ext,'\0',max_dim*max_thread_*sizeof(double));

                for ( int h_i = 0; h_i < nirrep_; h_i++){

                    int nmo_i = last_index_[h_i][0] - first_index_[h_i][0];
                    if ( nmo_i == 0 ) continue;

                    for ( int i_c = first_index_[h_i][0]; i_c < last_index_[h_i][0]; i_c++){

                        // gather all (:,ip) for this i

                        for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                            int ip_df  = df_ij_index(i_c,p_c);
                            int ip_scr = p_c - first_index_[h_p][p_class];

                            C_DCOPY(nQ_,&tei[ip_df],Qstride_,&scr_3index_1[ip_scr*nQ_],1);

                        } // end p_c loop

                        // gather all (:,iq) for this i

                        for ( int q_c = first_index_[h_p][q_class]; q_c < last_index_[h_p][q_class]; q_c++){

                            int iq_df  = df_ij_index(i_c,q_c);
                            int iq_scr = q_c - first_index_[h_p][q_class];

                            C_DCOPY(nQ_,&tei[iq_df],Qstride_,&scr_3index_2[iq_scr*nQ_],1);

                        } // end q_c loop

                        // calculate all (ip|iq) --> nmo_p x nmo_q

                       if ( nmo_p*nmo_q > 0 ) F_DGEMM('t','n',nmo_p,nmo_q,nQ_,1.0e0,&scr_3index_1[0],nQ_,&scr_3index_2[0],nQ_,1.0e0,&scr_4index_1[0],nmo_p);

                        if ( calc_q_diag[q_class*nirrep_ + h_p] ){

                            // calculate all (iq|iq)
                            if ( nmo_q > 0 ) F_DGEMM('t','n',nmo_q,nmo_q,nQ_,1.0e0,&scr_3index_2[0],nQ_,&scr_3index_2[0],nQ_,1.0e0,&scr_4index_2[0],nmo_q);

                        }

                        if ( calc_p_diag[p_class*nirrep_ + h_p] == 1 ){

                            for ( int p_i = 0; p_i < nmo_p; p_i++){

                                int thread_id = 0;

                                scr_diag_ext[thread_id*nmo_p+p_i] += C_DDOT(nQ_,&scr_3index_1[p_i*nQ_],1,&scr_3index_1[p_i*nQ_],1);                       

                            } // end p_i loop    

                        }

                    } // end i_c loop

                } // end h_i loop

                // update off-diagonal blocks of Fi (ca,ce,ae)

                for ( int q_c = first_index_[h_p][q_class]; q_c < last_index_[h_p][q_class]; q_c++){

                    int q_scr = q_c - first_index_[h_p][q_class];

                    for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                        int p_scr = p_c - first_index_[h_p][p_class];
                        int Fi_index = Lt_ij_oe_index(p_c,q_c,h_p);

                        Fi_[Fi_index] -= scr_4index_1[Full_ij_index(p_scr,q_scr,nmo_p)];

                    } // end p_c loop

                } // end p_c loop

                // update diagonal block (cc and aa)

                if ( calc_q_diag[q_class*nirrep_ + h_p] == 1 ){

                   for ( int q_c = first_index_[h_p][q_class]; q_c < last_index_[h_p][q_class]; q_c++){

                       int q_scr = q_c - first_index_[h_p][q_class];

                       for ( int p_c = q_c; p_c < last_index_[h_p][q_class]; p_c++){

                           int p_scr = p_c - first_index_[h_p][q_class];
                           int Fi_index = Lt_ij_oe_index(p_c,q_c,h_p);

                           Fi_[Fi_index] -= scr_4index_2[Full_ij_index(p_scr,q_scr,nmo_q)];

                       } // end p_c loop

                   } // end q_c loop

                   calc_q_diag[q_class*nirrep_ + h_p] = 0;

                }

                // update diagonal block (ee)

                if ( calc_p_diag[p_class*nirrep_ + h_p] ){

                   for ( int thread_reduce = 0; thread_reduce < max_thread_; thread_reduce++){

                       for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                           int p_scr = p_c - first_index_[h_p][p_class];
                           int Fi_index   = Lt_ij_oe_index(p_c,p_c,h_p);

                           Fi_[Fi_index] -= scr_diag_ext[thread_reduce*nmo_p+p_scr];

                       } // end p_c loop

                   } // end thread_reduce loop

                   calc_p_diag[p_class*nirrep_ + h_p] = 0;

                }

            } // end h_p loop

        } // end q_class loop

    } // end p_class loop

    free(scr_3index_1);
    free(scr_3index_2);
    free(scr_4index_1);
    free(scr_4index_2);
    free(scr_diag_ext);

}

void OrbitalOptimizer::Compute_Fi_J(double * tei, double * oei){

    // function to calculate Coulomb contribution to Fi
    // set up so that it can run in parallel (as long as thread_id is called with openmp function)
    // currently, only diagonal elements are calculate for e-e block (can be changed by the q_c_min variable)

    double * scr_3index_1 = (double *)malloc(nQ_*sizeof(double));
    memset((void*)scr_3index_1,'\0',nQ_*sizeof(double));

    double * scr_3index_2 = (double *)malloc(max_thread_*nQ_*sizeof(double));
    memset((void*)scr_3index_2,'\0',max_thread_*nQ_*sizeof(double));

    for ( int h_i = 0; h_i < nirrep_; h_i++){

        for ( int i_c = first_index_[h_i][0]; i_c < last_index_[h_i][0]; i_c++){

            int ii_df     = df_ij_index(i_c,i_c);
            int thread_id = 0;

            C_DAXPY(nQ_,2.0e0,&tei[ii_df],Qstride_,&scr_3index_2[thread_id*nQ_],1);

        }

    }

    for (int thread_reduce = 0; thread_reduce < max_thread_; thread_reduce++){

         C_DAXPY(nQ_,1.0e0,&scr_3index_2[thread_reduce*nQ_],1,&scr_3index_1[0],1);

    }

    for ( int p_class = 0; p_class < 3; p_class++ ){

        for ( int q_class = 0; q_class <= p_class; q_class++){

            for ( int h_p = 0; h_p < nirrep_; h_p++){

                for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                    int q_c_min = first_index_[h_p][q_class];
                    if ( q_class == 2 ) q_c_min = p_c;

                    int q_c_max = last_index_[h_p][q_class];
                    if ( p_class == q_class ) q_c_max = p_c + 1;

                    for ( int q_c = q_c_min; q_c < q_c_max; q_c++){

                        int pq_df = df_ij_index(p_c,q_c);

                        int oei_index = Full_GindMap_c_[Full_ij_index(p_c,q_c,nmo_)];
                        int Fi_index  = Lt_ij_oe_index(p_c,q_c,h_p);

                        Fi_[Fi_index] = oei[oei_index] + C_DDOT(nQ_,&tei[pq_df],Qstride_,&scr_3index_1[0],1);

                    } // end q_c loop 

                } // end p_c loop

            } // end h_p loop

        } // end q_class loop

    }  // end p_class loop


    free(scr_3index_1);
    free(scr_3index_2);

}

void OrbitalOptimizer::Initialize_scratch_G(){
    int Q_dim = 0;
    for ( int h = 0; h< nirrep_; h++) Q_dim += amopi_[h] * nmopi_[h];

    int fock_dim = 0;
    for ( int h = 0; h < nirrep_; h++) fock_dim += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;

    memset((void*)Fi_,'\0',fock_dim*sizeof(double));

    memset((void*)Fa_,'\0',fock_dim*sizeof(double));

    memset((void*)Q_,'\0',Q_dim*sizeof(double));

    memset((void*)Z_,'\0',Q_dim*sizeof(double));

}

void OrbitalOptimizer::Alloc_Gradient(){

    int fock_dim = 0;
    for ( int h = 0; h < nirrep_; h++) fock_dim += nmopi_[h] * ( nmopi_[h] + 1 ) / 2;

    Fi_ = (double *)malloc(fock_dim*sizeof(double));
    memset((void*)Fi_,'\0',fock_dim*sizeof(double));

    Fa_ = (double *)malloc(fock_dim*sizeof(double));
    memset((void*)Fa_,'\0',fock_dim*sizeof(double));

    int Q_dim = 0;
    for ( int h = 0; h< nirrep_; h++) Q_dim += amopi_[h] * nmopi_[h];

    Q_ = (double *)malloc(Q_dim*sizeof(double));
    memset((void*)Q_,'\0',Q_dim*sizeof(double));

    Z_ = (double *)malloc(Q_dim*sizeof(double));
    memset((void*)Z_,'\0',Q_dim*sizeof(double));

    Grad_IndMap_c_ = (int *)malloc(nmo_*(nmo_+1)/2*sizeof(int));
    memset((void*)Grad_IndMap_c_,'\0',nmo_*(nmo_+1)/2*sizeof(int));

    Nrot_ = 0;

    for ( int p_class = 0; p_class < 3; p_class++ ){

        int q_class_max = p_class;
        if ( active_active_rotations_ && p_class == 1 ) q_class_max = p_class + 1;

        for ( int q_class = 0; q_class < q_class_max; q_class++ ){

            for ( int h = 0; h < nirrep_; h++ ){

                int nmo_p = last_index_[h][p_class] - first_index_[h][p_class];

                int nmo_q = last_index_[h][q_class] - first_index_[h][q_class];

                if ( p_class == q_class ){

                   Nrot_ += nmo_p * ( nmo_p - 1 ) / 2;

                }

                else {

                   Nrot_ += nmo_p * nmo_q;

                }

            } // end h loop

        } // emd q_class loop       

    } // end p_class loop

    orbital_gradient_ = (double *)malloc(Nrot_*sizeof(double));
    memset((void*)orbital_gradient_,'\0',Nrot_*sizeof(double));

}

void OrbitalOptimizer::Dealloc_Gradient(){

     free(Fi_);
     free(Fa_);
     free(Q_);
     free(Z_);
     free(Grad_IndMap_c_);
     free(orbital_gradient_);

}

void OrbitalOptimizer::Printfock(double * fock){

    for ( int h=0; h < nirrep_; h++){

        for ( int p_class = 0; p_class < 3; p_class++){

            for ( int q_class=0; q_class < p_class+1; q_class++){

                for ( int p_c=first_index_[h][p_class]; p_c<last_index_[h][p_class]; p_c++){

                    for ( int q_c=first_index_[h][q_class]; q_c<last_index_[h][q_class]; q_c++){

                        if ( q_c > p_c ) continue;
                        if ( p_class == 2 && q_class == 2 && p_c != q_c ) continue;
                        int Fi_index = Lt_ij_oe_index(p_c,q_c,h);

                        if ( fock[Fi_index] == 0.0e0 ) continue;
                        printf("%1i    %1i %1i     %3i %3i    %20.12e\n",h,p_class,q_class,p_c,q_c,fock[Fi_index]);

                    }

                }

            }

        }

    }

}



}
