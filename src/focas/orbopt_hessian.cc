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

void OrbitalOptimizer::ComputeH_diag(double * d2, double * d1, double * tei){

    memset((void*)diagonal_orbital_hessian_,'\0',Nrot_*sizeof(double));

    ComputeH_diag_ad(d2, d1, tei);

    ComputeH_diag_ed(tei);

    if ( active_active_rotations_ ) ComputeH_diag_aa(d2, d1, tei);

    ComputeH_diag_ea(d2, d1, tei);

}

void OrbitalOptimizer::CleanH_diag(){

    for ( int pq = 0; pq < Nrot_; pq++){

        double Hval = diagonal_orbital_hessian_[pq];

        if ( Hval < 0.0e0   ) Hval = - Hval;
        if ( Hval <= 1.0e-2 ) Hval = 1.0e-2;

        diagonal_orbital_hessian_[pq] = Hval;

    }

}

void OrbitalOptimizer::PrintH_diag(){

    printf("\n");
    printf("Diagonal elements of H (a-d)\n");
    printf("\n");
    for ( int h = 0; h < nirrep_; h++){

        for ( int p = first_index_[h][1]; p < last_index_[h][1]; p++){

            for ( int q = first_index_[h][0]; q < last_index_[h][0]; q++){

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(p,q)];
                printf("%i     %i %i     %20.13e\n",h,p,q,diagonal_orbital_hessian_[grad_ind]);

            }

        }

    }

    printf("\n");
    printf("Diagonal elements of H (e-d)\n");
    printf("\n");
    for ( int h = 0; h < nirrep_; h++){

        for ( int p = first_index_[h][2]; p < last_index_[h][2]; p++){

            for ( int q = first_index_[h][0]; q < last_index_[h][0]; q++){

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(p,q)];
                printf("%i     %i %i     %20.13e\n",h,p,q,diagonal_orbital_hessian_[grad_ind]);

            }

        }

    }

    if ( active_active_rotations_ ){

        printf("\n");
        printf("Diagonal elements of H (a-a)\n");
        printf("\n");
        for ( int h = 0; h < nirrep_; h++){

            for ( int p = first_index_[h][1]; p < last_index_[h][1]; p++){

                for ( int q = first_index_[h][1]; q < p; q++){

                    int grad_ind = Grad_IndMap_c_[Lt_ij_index(p,q)];
                    printf("%i     %i %i     %20.13e\n",h,p,q,diagonal_orbital_hessian_[grad_ind]);

                }

            }

        }


    }

    printf("\n");
    printf("Diagonal elements of H (e-a)\n");
    printf("\n");
    for ( int h = 0; h < nirrep_; h++){

        for ( int p = first_index_[h][2]; p < last_index_[h][2]; p++){

            for ( int q = first_index_[h][1]; q < last_index_[h][1]; q++){

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(p,q)];
                printf("%i     %i %i     %20.13e\n",h,p,q,diagonal_orbital_hessian_[grad_ind]);

            }

        }

    }

}

void OrbitalOptimizer::ComputeH_diag_ad(double * d2, double * d1, double * tei){

    // Diagonal Hessian elements
    // exact -- Eq. 4.6c && approximate Eq. 4.7c in Gordon, TCA, 97, 88-95 (1997) 

    for ( int h = 0; h < nirrep_; h++ ){

        if ( amopi_[h] * rstcpi_[h] == 0 ) continue;

        for ( int t_c = first_index_[h][1]; t_c < last_index_[h][1]; t_c++ ){

            for ( int i_c = first_index_[h][0]; i_c < last_index_[h][0]; i_c++ ){

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(t_c,i_c)];
                int tt_den   = Active_GindMap_c_[Full_ij_index(t_c,t_c,nmo_)];

                diagonal_orbital_hessian_[grad_ind] =
                          
                    2.0e0 * (

                        2.0e0 * (

                              Fi_[ Lt_ij_oe_index(t_c,t_c,h) ]
                            + Fa_[ Lt_ij_oe_index(t_c,t_c,h) ]
                            - Fi_[ Lt_ij_oe_index(i_c,i_c,h) ]
                            - Fa_[ Lt_ij_oe_index(i_c,i_c,h) ]

                        )

                        + d1[ tt_den ] * Fi_[ Lt_ij_oe_index(i_c,i_c,h) ]
                        - Q_[ Q_index(t_c,t_c,h) ]
                        - Z_[ Q_index(t_c,t_c,h) ]

                    );

                if ( !exact_diagonal_hessian_ ) {

                    diagonal_orbital_hessian_[grad_ind] += 2.0e0 * d1[ tt_den ] * Fa_[ Lt_ij_oe_index(i_c,i_c,h) ];

                }

            } // end i_c loop

        } // end t_c loop

        if ( exact_diagonal_hessian_ ) {

            tei_terms_H_ad(d2, d1, tei, h);

        }

    } // end h loop

}

void OrbitalOptimizer::tei_terms_H_ad(double * d2, double * d1, double * tei, int h){

    // allocate scratch arrays

    int ngem_ts    = act_gempi_[0];
    int max_doc    = Max_value(rstcpi_,nirrep_);
    int max_4index = ngem_ts; if ( max_doc > ngem_ts ) max_4index = max_doc;

    double * scr_4index_1 = (double *)malloc(max_4index*max_4index*sizeof(double));
    double * scr_3index_1 = (double *)malloc(ngem_ts*nQ_*sizeof(double));
    double * scr_3index_2 = (double *)malloc(max_doc*nQ_*sizeof(double));


    // *** Coulomb terms ***

    // gather all (ii|:) for this h

    int ngem_ii = rstcpi_[h];

    for ( int i_c = first_index_[h][0]; i_c < last_index_[h][0]; i_c++ ){

        int ii_scr = ( i_c - first_index_[h][0] ) * nQ_;
        int ii_df  = df_ij_index(i_c,i_c);

        C_DCOPY(nQ_,&tei[ii_df],Qstride_,&scr_3index_2[ii_scr],1);

    }

    // gather all (uv|:) for all h_u=h_v

    for ( int h_u = 0; h_u < nirrep_; h_u++){

        // gather all (tu|:) for this h_u

        for ( int u_c = first_index_[h_u][1]; u_c < last_index_[h_u][1]; u_c++ ){

            for ( int v_c = first_index_[h_u][1]; v_c < u_c + 1; v_c++ ){

                int uv_df  = df_ij_index(u_c,v_c);
                int uv_scr = Active_GindMap_c_[Full_ij_index(u_c,v_c,nmo_)] * nQ_;

                C_DCOPY(nQ_,&tei[uv_df],Qstride_,&scr_3index_1[uv_scr],1);

            }

        }

    }

    // Calculate all (uv|ii) for all h_u and h

   if ( ngem_ts > 0 ) F_DGEMM('t','n',ngem_ts,ngem_ii,nQ_,1.0e0,&scr_3index_1[0],nQ_,&scr_3index_2[0],nQ_,0.0e0,&scr_4index_1[0],ngem_ts);

    // *** update diagonal Hessian

    for ( int h_u = 0; h_u < nirrep_; h_u++){

        // 2*d2(tt|uv)*tei(uv|ii) for all u&v

        for ( int t_c = first_index_[h][1]; t_c < last_index_[h][1]; t_c++){

            int tt_den = Active_GindMap_c_[Full_ij_index(t_c,t_c,nmo_)];

            for ( int i_c = first_index_[h][0]; i_c < last_index_[h][0]; i_c++){

                int ii = i_c - first_index_[h][0];

                double Hval = 0.0e0;

                // 4*d2(tt|uv)*tei(uv|ii) since only LT elements are accessed

                for ( int u_c = first_index_[h_u][1]; u_c < last_index_[h_u][1]; u_c++){

                    for ( int v_c = first_index_[h_u][1]; v_c < u_c; v_c++){

                        int uv_den   = Active_GindMap_c_[Full_ij_index(u_c,v_c,nmo_)];
                        int ttuv_den = Lt_ij_index(tt_den,uv_den);
                        int iiuv_int = Full_ij_index(uv_den,ii,ngem_ts);

                        Hval += 4.0e0 * d2[ttuv_den] * scr_4index_1[iiuv_int];

                    }

                    // 2*d2(tt|uu)*tei(uu|ii)      

                    int uu_den   = Active_GindMap_c_[Full_ij_index(u_c,u_c,nmo_)];
                    int ttuu_den = Lt_ij_index(tt_den,uu_den);
                    int iiuu_int = Full_ij_index(uu_den,ii,ngem_ts);

                    Hval += 2.0e0 * d2[ttuu_den] * scr_4index_1[iiuu_int];

                }

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(t_c,i_c)];

                diagonal_orbital_hessian_[grad_ind] += Hval;

            }

        }

        if ( h_u != h ) continue;

        // 4*tei(tu|ii) * ( d1(t,u) - delta(t,u) )

        for ( int t_c = first_index_[h][1]; t_c < last_index_[h][1]; t_c++){

            for ( int i_c = first_index_[h][0]; i_c < last_index_[h][0]; i_c++){

                int ii = i_c - first_index_[h][0];


                int tt_den      = Active_GindMap_c_[Full_ij_index(t_c,t_c,nmo_)];
                int iitt_int    = Full_ij_index(tt_den,ii,ngem_ts);
                double Hval     = -scr_4index_1[iitt_int];

                for ( int u_c = first_index_[h][1]; u_c < last_index_[h][1]; u_c++){

                    int tu_den   = Active_GindMap_c_[Full_ij_index(t_c,u_c,nmo_)];
                    int iitu_int = Full_ij_index(tu_den,ii,ngem_ts);

                    Hval += d1[tu_den] * scr_4index_1[iitu_int];

                }

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(t_c,i_c)];

                diagonal_orbital_hessian_[grad_ind] += 4.0e0 * Hval;

            }

        }

    }

    // exchange terms

    for ( int i_c = first_index_[h][0]; i_c < last_index_[h][0]; i_c++ ){

        // gather all (ui|:) for all u and this i

        for ( int h_u = 0; h_u < nirrep_; h_u++ ){

            for ( int u_c = first_index_[h_u][1]; u_c < last_index_[h_u][1]; u_c++){

                int ui_df  = df_ij_index(u_c,i_c);
                int ui_scr = ( u_c - nrstc_ ) * nQ_;

                C_DCOPY(nQ_,&tei[ui_df],Qstride_,&scr_3index_1[ui_scr],1);

            }

        }

/// GG note to self ... I should rewrite this so that only required integrals are calculated here (see exchange terms in ea function)

        // Calculate all (ui|vi) for all u, v and this i

        if ( amo_ > 0 ) F_DGEMM('t','n',amo_,amo_,nQ_,1.0e0,&scr_3index_1[0],nQ_,&scr_3index_1[0],nQ_,0.0e0,&scr_4index_1[0],amo_);

        // update Hessian elements

        for ( int h_u = 0; h_u < nirrep_; h_u++){

            int d2_offset = d2_irrep_offset_[GemSym(h,h_u)];

            // 4 * d2(tv|tu) * tei(ui|vi)            

            for ( int t_c = first_index_[h][1]; t_c < last_index_[h][1]; t_c++){

                double Hval = 0.0e0;

                for ( int u_c = first_index_[h_u][1]; u_c < last_index_[h_u][1]; u_c++){

                    for ( int v_c = first_index_[h_u][1]; v_c < last_index_[h_u][1]; v_c++){

                        int tu_den   = Active_GindMap_c_[Full_ij_index(t_c,u_c,nmo_)];
                        int tv_den   = Active_GindMap_c_[Full_ij_index(t_c,v_c,nmo_)];
                        int tutv_den = Lt_ij_index(tu_den,tv_den) + d2_offset;

                        int ui_scr   = u_c - nrstc_;
                        int vi_scr   = v_c - nrstc_;
                        int uivi_scr = Full_ij_index(ui_scr,vi_scr,amo_);

                        Hval += 4.0e0 * d2[tutv_den] * scr_4index_1[uivi_scr];

                    }

                }

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(t_c,i_c)];

                diagonal_orbital_hessian_[grad_ind] += Hval;

            }

            if ( h != h_u ) continue;

            // 12 * tei(ui|ti) * ( delta(t,u) - d1(t,u) )

            for ( int t_c = first_index_[h][1]; t_c < last_index_[h][1]; t_c++){

                int ti_scr   = t_c - nrstc_;
                int titi_scr = Full_ij_index(ti_scr,ti_scr,amo_);
                double Hval  = scr_4index_1[titi_scr];

                for ( int u_c = first_index_[h][1]; u_c < last_index_[h][1]; u_c++){

                    int ui_scr   = u_c - nrstc_;
                    int ti_scr   = t_c - nrstc_;
                    int uiti_scr = Full_ij_index(ui_scr,ti_scr,amo_);

                    int tu_den   = Active_GindMap_c_[Full_ij_index(t_c,u_c,nmo_)];

                    Hval        -= d1[tu_den] * scr_4index_1[uiti_scr];

                }

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(t_c,i_c)];

                diagonal_orbital_hessian_[grad_ind] += 12.0e0 * Hval;

            }

        }

    }

    // deallocate scratch arrays

    free(scr_3index_1);
    free(scr_3index_2);
    free(scr_4index_1);

}

void OrbitalOptimizer::ComputeH_diag_ed(double * tei){

    // Diagonal Hessian elements
    // exact -- Eq. 4.6a && approximate Eq. 4.7a in Gordon, TCA, 97, 88-95 (1997) 

    for ( int h = 0; h < nirrep_; h++ ){

        for ( int a_c = first_index_[h][2]; a_c < last_index_[h][2]; a_c++ ){

            for ( int i_c = first_index_[h][0]; i_c < last_index_[h][0]; i_c++ ){

                double Hval = Fi_[ Lt_ij_oe_index(a_c,a_c,h) ]
                            + Fa_[ Lt_ij_oe_index(a_c,a_c,h) ]
                            - Fi_[ Lt_ij_oe_index(i_c,i_c,h) ]
                            - Fa_[ Lt_ij_oe_index(i_c,i_c,h) ];

                if ( exact_diagonal_hessian_ ) {

                    int ia_df = df_ij_index(i_c,a_c);
                    int ii_df = df_ij_index(i_c,i_c);
                    int aa_df = df_ij_index(a_c,a_c);

                    Hval += 3.0e0 * C_DDOT(nQ_,&tei[ia_df],Qstride_,&tei[ia_df],Qstride_)
                                  - C_DDOT(nQ_,&tei[ii_df],Qstride_,&tei[aa_df],Qstride_);

                }

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(a_c,i_c)];

                diagonal_orbital_hessian_[grad_ind] = 4.0e0 * Hval;

            } // end i_c loop

        } // end t_c loop

    } // end h loop

}

void OrbitalOptimizer::ComputeH_diag_aa(double * d2, double * d1, double * tei){

    /*
    2-e terms adapted from Eq. B27 in Shepard, 

    H(tu|tu) += sum_rs [ 2 * g(tr|ts) * d2(ur|us) 
                       + 2 * g(ur|us) * d2(tr|ts) 
                       - 4 * g(ur|ts) * d2(ur|ts) 
                       +     g(tt|rs) * d2(uu|rs)  
                       +     g(uu|rs) * d2(tt|rs)    
                       - 2 * g(tu|rs) * d2(tu|rs) 
    */

    // allocate scratch arrays for 

    int max_act = Max_value(amopi_,nirrep_);
    int tmp     = 0;
    int h_max   = 0;

    for ( int h = 0; h < nirrep_; h++) if ( amopi_[h] == max_act ) h_max = h;

    for ( int h = 0; h < nirrep_; h++){

        if ( h == h_max ) continue;
        if ( tmp < amopi_[h] ) tmp = amopi_[h];

    }

    int dim_scr_2 = max_act * ( max_act + 1 ) / 2;
    int dim_scr_1 = max_act *  tmp;
    if ( dim_scr_1 < dim_scr_2 ) dim_scr_1 = dim_scr_2;

    double * scr_4index_1 = (double *)malloc(dim_scr_1*dim_scr_1*sizeof(double));

    // storage for (rs|:)
    double * scr_3index_1 = (double *)malloc(dim_scr_1*nQ_*sizeof(double));

    // storage for (tr|:)
    double * scr_3index_2 = (double *)malloc(dim_scr_2*nQ_*sizeof(double));

    // 1-e terms

    for ( int h = 0; h < nirrep_; h++ ){

        for ( int t_c = first_index_[h][1]; t_c < last_index_[h][1]; t_c++ ){

            for ( int u_c = first_index_[h][1]; u_c < t_c; u_c++ ){

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(t_c,u_c)];
                int tu_den   = Active_GindMap_c_[Full_ij_index(t_c,u_c,nmo_)];
                int uu_den   = Active_GindMap_c_[Full_ij_index(u_c,u_c,nmo_)];
                int tt_den   = Active_GindMap_c_[Full_ij_index(t_c,t_c,nmo_)];


                diagonal_orbital_hessian_[grad_ind] =
                                  2.0e0 * (
                                        - ( Q_[Q_index(t_c,t_c,h)] \
                                          + Z_[Q_index(t_c,t_c,h)] \
                                          + Q_[Q_index(u_c,u_c,h)] \
                                          + Z_[Q_index(u_c,u_c,h)] \
                                          + 2.0e0 * d1[tu_den] * Fi_[ Lt_ij_oe_index(t_c,u_c,h) ] ) \
                                        + d1[tt_den] * Fi_[ Lt_ij_oe_index(u_c,u_c,h) ] \
                                        + d1[uu_den] * Fi_[ Lt_ij_oe_index(t_c,t_c,h) ] \
                                          );

            } // end u_c loop 

        } // end t_c loop

    } // end h loop

    // 2-e coulomb terms

    for ( int h_r = 0; h_r < nirrep_; h_r++ ){

        int ngem_rs = amopi_[h_r] * ( amopi_[h_r] + 1 ) / 2;
        if ( ngem_rs == 0 ) continue;

        // gather all (rs|:) for this h_r

        for ( int r_c = first_index_[h_r][1]; r_c < last_index_[h_r][1]; r_c++ ){

            for ( int s_c = first_index_[h_r][1]; s_c < r_c + 1; s_c++ ){

                int r_i    = r_c - first_index_[h_r][1];
                int s_i    = s_c - first_index_[h_r][1];
                int rs_scr = Lt_ij_index(r_i,s_i) * nQ_;
                int rs_df  = df_ij_index(r_c,s_c);

                C_DCOPY(nQ_,&tei[rs_df],Qstride_,&scr_3index_1[rs_scr],1);

            }

        }

        for ( int h_t = 0; h_t < nirrep_; h_t++){

            int ngem_tu = amopi_[h_t] * ( amopi_[h_t] + 1 ) / 2;

            if ( ngem_tu != 0 &&  ngem_rs != 0){

                // *************************
                // *** 2-e coulomb terms ***
                // *************************

                // gather all (tu|:) for this h_t != h_r

                if ( h_r != h_t ) {

                    for ( int t_c = first_index_[h_t][1]; t_c < last_index_[h_t][1]; t_c++ ){

                        for ( int u_c = first_index_[h_t][1]; u_c < t_c + 1; u_c++ ){

                            int t_i    = t_c - first_index_[h_t][1];
                            int u_i    = u_c - first_index_[h_t][1];
                            int tu_scr = Lt_ij_index(t_i,u_i) * nQ_;
                            int tu_df  = df_ij_index(t_c,u_c);

                            C_DCOPY(nQ_,&tei[tu_df],Qstride_,&scr_3index_2[tu_scr],1);

                        }

                    }

                }

                // compute all (rs|tu) for this h_t and h_r

                if ( h_r != h_t ) {

                    F_DGEMM('t','n',ngem_rs,ngem_tu,nQ_,1.0e0,&scr_3index_1[0],nQ_,\
                               &scr_3index_2[0],nQ_,0.0e0,&scr_4index_1[0],ngem_rs);

                }

                else {

                    F_DGEMM('t','n',ngem_rs,ngem_tu,nQ_,1.0e0,&scr_3index_1[0],nQ_,\
                               &scr_3index_1[0],nQ_,0.0e0,&scr_4index_1[0],ngem_rs);

                }

                // update Hessian elements
                /*
                    H(tu|tu) += sum_rs [ g(tt|rs) * d2(uu|rs)  
                                       + g(uu|rs) * d2(tt|rs)    
                                   - 2 * g(tu|rs) * d2(tu|rs)
                */

                for ( int t_c = first_index_[h_t][1]; t_c < last_index_[h_t][1]; t_c++ ){

                    for ( int u_c = first_index_[h_t][1]; u_c < t_c ; u_c++ ){

                        double Hval = 0.0e0;

                        for ( int r_c = first_index_[h_r][1]; r_c < last_index_[h_r][1]; r_c++ ){

                            for ( int s_c = first_index_[h_r][1]; s_c < last_index_[h_r][1]; s_c++ ){

                                int t_i      = t_c - first_index_[h_t][1];
                                int u_i      = u_c - first_index_[h_t][1];
                                int r_i      = r_c - first_index_[h_r][1];
                                int s_i      = s_c - first_index_[h_r][1];

                                int tt_scr   = Lt_ij_index(t_i,t_i);
                                int uu_scr   = Lt_ij_index(u_i,u_i);
                                int tu_scr   = Lt_ij_index(t_i,u_i);
                                int rs_scr   = Lt_ij_index(r_i,s_i);

                                int tt_den   = Active_GindMap_c_[Full_ij_index(t_c,t_c,nmo_)];
                                int uu_den   = Active_GindMap_c_[Full_ij_index(u_c,u_c,nmo_)];
                                int tu_den   = Active_GindMap_c_[Full_ij_index(t_c,u_c,nmo_)];
                                int rs_den   = Active_GindMap_c_[Full_ij_index(r_c,s_c,nmo_)];

                                int ttrs_int = Full_ij_index(rs_scr,tt_scr,ngem_rs);
                                int uurs_int = Full_ij_index(rs_scr,uu_scr,ngem_rs);
                                int turs_int = Full_ij_index(rs_scr,tu_scr,ngem_rs);

                                int ttrs_den = Lt_ij_index(tt_den,rs_den);
                                int uurs_den = Lt_ij_index(uu_den,rs_den);
                                int turs_den = Lt_ij_index(tu_den,rs_den);

                                Hval +=         tei[ttrs_int] * d2[uurs_den] \
                                              + tei[uurs_int] * d2[ttrs_den] \
                                      - 2.0e0 * tei[turs_int] * d2[turs_den];

                            }

                        }

                        int grad_ind = Grad_IndMap_c_[Lt_ij_index(t_c,u_c)];

                        diagonal_orbital_hessian_[grad_ind] += Hval;

                    }

                }

            } // end ngem_tu & ngem_rs (coulomb terms)

            // **************************
            // *** 2-e exchange terms ***
            // **************************

            int nmo_h_t  = amopi_[h_t];
            int nmo_h_r  = amopi_[h_r];
            int h_tr     = GemSym(h_t,h_r);
            int d2_h_off = d2_irrep_offset_[h_tr];

            if ( nmo_h_t * nmo_h_r != 0 ){

                if ( h_r != h_t ) {

                    int ngem_tr = nmo_h_t*nmo_h_r;

                    // collect all (rt|:) integrals for this h_t & h_r

                    for ( int t_c = first_index_[h_t][1]; t_c < last_index_[h_t][1]; t_c++ ){

                        for ( int r_c = first_index_[h_r][1]; r_c < last_index_[h_r][1]; r_c++ ){

                            int t_i    = t_c - first_index_[h_t][1];
                            int r_i    = r_c - first_index_[h_r][1];
                            int tr_scr = Full_ij_index(r_i,t_i,nmo_h_r)  * nQ_;
                            int tr_df  = df_ij_index(t_c,r_c);

                            C_DCOPY(nQ_,&tei[tr_df],Qstride_,&scr_3index_1[tr_scr],1);

                        }

                    }

                    // calculate all integrals (tr|us)

                    F_DGEMM('t','n',ngem_tr,ngem_tr,nQ_,1.0e0,&scr_3index_1[0],nQ_,\
                               &scr_3index_1[0],nQ_,0.0e0,&scr_4index_1[0],ngem_tr);

                    // update Hessian matrix elements

                    /*
                    H(tu|tu) += sum_rs [ 2 * g(tr|ts) * d2(ur|us) 
                                       + 2 * g(ur|us) * d2(tr|ts) 
                                       - 4 * g(ur|ts) * d2(ur|ts)
                    */

                    for ( int t_c = first_index_[h_t][1]; t_c < last_index_[h_t][1]; t_c++ ){

                        for ( int u_c = first_index_[h_t][1]; u_c < t_c ; u_c++ ){

                            double Hval = 0.0e0;

                            for ( int r_c = first_index_[h_r][1]; r_c < last_index_[h_r][1]; r_c++ ){

                                for ( int s_c = first_index_[h_r][1]; s_c < last_index_[h_r][1]; s_c++ ){

                                    int t_i      = t_c - first_index_[h_t][1];
                                    int u_i      = u_c - first_index_[h_t][1];
                                    int r_i      = r_c - first_index_[h_r][1];
                                    int s_i      = s_c - first_index_[h_r][1];

                                    int tr_scr   = Full_ij_index(r_i,t_i,nmo_h_r);
                                    int ts_scr   = Full_ij_index(s_i,t_i,nmo_h_r);
                                    int ur_scr   = Full_ij_index(r_i,u_i,nmo_h_r);
                                    int us_scr   = Full_ij_index(s_i,u_i,nmo_h_r);

                                    int tr_den   = Active_GindMap_c_[Full_ij_index(t_c,r_c,nmo_)];
                                    int ts_den   = Active_GindMap_c_[Full_ij_index(t_c,s_c,nmo_)];
                                    int ur_den   = Active_GindMap_c_[Full_ij_index(u_c,r_c,nmo_)];
                                    int us_den   = Active_GindMap_c_[Full_ij_index(u_c,s_c,nmo_)];

                                    int trts_int = Full_ij_index(tr_scr,ts_scr,ngem_tr);
                                    int urus_int = Full_ij_index(ur_scr,us_scr,ngem_tr);
                                    int urts_int = Full_ij_index(ur_scr,ts_scr,ngem_tr);

                                    int trts_den = Lt_ij_index(tr_den,ts_den) + d2_h_off;
                                    int urus_den = Lt_ij_index(ur_den,us_den) + d2_h_off;
                                    int urts_den = Lt_ij_index(ur_den,ts_den) + d2_h_off;

                                    Hval +=         tei[trts_int] * d2[urus_den] \
                                                  + tei[urus_int] * d2[trts_den] \
                                          - 2.0e0 * tei[urts_int] * d2[urts_den];

                                }

                            }

                            int grad_ind = Grad_IndMap_c_[Lt_ij_index(t_c,u_c)];

                            diagonal_orbital_hessian_[grad_ind] += 2.0e0 * Hval;

                        }

                    }

                }

                else {

                    int ngem_tr = nmo_h_t * ( nmo_h_t + 1 ) / 2;

                    // (tr|us) integrals are already stored scr_4index_1

                    for ( int t_c = first_index_[h_t][1]; t_c < last_index_[h_t][1]; t_c++ ){

                        for ( int u_c = first_index_[h_t][1]; u_c < t_c ; u_c++ ){

                            double Hval = 0.0e0;

                            for ( int r_c = first_index_[h_r][1]; r_c < last_index_[h_r][1]; r_c++ ){

                                for ( int s_c = first_index_[h_r][1]; s_c < last_index_[h_r][1]; s_c++ ){

                                    int t_i      = t_c - first_index_[h_t][1];
                                    int u_i      = u_c - first_index_[h_t][1];
                                    int r_i      = r_c - first_index_[h_r][1];
                                    int s_i      = s_c - first_index_[h_r][1];

                                    int tr_scr   = Lt_ij_index(t_i,r_i);
                                    int ts_scr   = Lt_ij_index(t_i,s_i);
                                    int ur_scr   = Lt_ij_index(u_i,r_i);
                                    int us_scr   = Lt_ij_index(u_i,s_i);

                                    int tr_den   = Active_GindMap_c_[Full_ij_index(t_c,r_c,nmo_)];
                                    int ts_den   = Active_GindMap_c_[Full_ij_index(t_c,s_c,nmo_)];
                                    int ur_den   = Active_GindMap_c_[Full_ij_index(u_c,r_c,nmo_)];
                                    int us_den   = Active_GindMap_c_[Full_ij_index(u_c,s_c,nmo_)];

                                    int trts_int = Full_ij_index(tr_scr,ts_scr,ngem_tr);
                                    int urus_int = Full_ij_index(ur_scr,us_scr,ngem_tr);
                                    int urts_int = Full_ij_index(ur_scr,ts_scr,ngem_tr);

                                    int trts_den = Lt_ij_index(tr_den,ts_den);
                                    int urus_den = Lt_ij_index(ur_den,us_den);
                                    int urts_den = Lt_ij_index(ur_den,ts_den);

                                    Hval +=         tei[trts_int] * d2[urus_den] \
                                                  + tei[urus_int] * d2[trts_den] \
                                          - 2.0e0 * tei[urts_int] * d2[urts_den];

                                }

                            }

                            int grad_ind = Grad_IndMap_c_[Lt_ij_index(t_c,u_c)];

                            diagonal_orbital_hessian_[grad_ind] += 2.0e0 * Hval;

                        }

                    }

                }

            } // end ngem_tr (exchange terms)

        } // end h_t

    } // end h_r

    free(scr_3index_1);
    free(scr_3index_2);
    free(scr_4index_1);

}

void OrbitalOptimizer::ComputeH_diag_ea(double * d2, double * d1, double * tei){

    // Diagonal Hessian elements
    // exact -- Eq. 4.6b && approximate Eq. 4.7b in Gordon, TCA, 97, 88-95 (1997) 

    for ( int h = 0; h < nirrep_; h++ ){

        if ( amopi_[h] * rstvpi_[h] == 0 ) continue;

        for ( int a_c = first_index_[h][2]; a_c < last_index_[h][2]; a_c++ ){

            for ( int t_c = first_index_[h][1]; t_c < last_index_[h][1]; t_c++ ){

                int tt_den   = Active_GindMap_c_[Full_ij_index(t_c,t_c,nmo_)];
                int grad_ind = Grad_IndMap_c_[Lt_ij_index(a_c,t_c)];

                diagonal_orbital_hessian_[grad_ind] =

                    2.0e0 * (

                              d1[ tt_den ] * Fi_[ Lt_ij_oe_index(a_c,a_c,h) ]
                            - Q_[ Q_index(t_c,t_c,h) ]
                            - Z_[ Q_index(t_c,t_c,h) ]

                            );

                if ( exact_diagonal_hessian_ ) continue;

                diagonal_orbital_hessian_[grad_ind] += 2.0e0 * d1[ tt_den ] * Fa_[ Lt_ij_oe_index(a_c,a_c,h) ];

            } // end i_c loop

        } // end t_c loop

        if ( !exact_diagonal_hessian_ ) continue;

        tei_terms_H_ea(d2, tei, h);

    } // end h loop

}

void OrbitalOptimizer::tei_terms_H_ea(double * d2, double * tei, int h){

    // allocate scratch arrays

    int ngem_ts    = act_gempi_[0];
    int max_ext    = Max_value(rstvpi_,nirrep_);
    int max_4index = ngem_ts; if ( max_ext > ngem_ts ) max_4index = max_ext;

    double * scr_4index_1 = (double *)malloc(max_4index*max_4index*sizeof(double));
    double * scr_3index_1 = (double *)malloc(ngem_ts*nQ_*sizeof(double));
    double * scr_3index_2 = (double *)malloc(max_ext*nQ_*sizeof(double));

    // Coulomb terms

    int ngem_aa = rstvpi_[h];

    for ( int a_c = first_index_[h][2]; a_c < last_index_[h][2]; a_c++ ){

        int aa_scr = ( a_c - first_index_[h][2] ) * nQ_;
        int aa_df  = df_ij_index(a_c,a_c);

        C_DCOPY(nQ_,&tei[aa_df],Qstride_,&scr_3index_2[aa_scr],1);

    }

    // gather all (uv|:) for all h_u=h_v

    for ( int h_u = 0; h_u < nirrep_; h_u++){

        // gather all (tu|:) for this h_u

        for ( int u_c = first_index_[h_u][1]; u_c < last_index_[h_u][1]; u_c++ ){

            for ( int v_c = first_index_[h_u][1]; v_c < u_c + 1; v_c++ ){

                int uv_df  = df_ij_index(u_c,v_c);
                int uv_scr = Active_GindMap_c_[Full_ij_index(u_c,v_c,nmo_)] * nQ_;

                C_DCOPY(nQ_,&tei[uv_df],Qstride_,&scr_3index_1[uv_scr],1);

            }

        }

    }

    // Calculate all (uv|aa) for all h_u and h

    if ( ngem_ts > 0 ) F_DGEMM('t','n',ngem_ts,ngem_aa,nQ_,1.0e0,&scr_3index_1[0],nQ_,&scr_3index_2[0],nQ_,0.0e0,&scr_4index_1[0],ngem_ts);

    // *** update diagonal Hessian

    for ( int h_u = 0; h_u < nirrep_; h_u++){

        // 2*d2(tt|uv)*tei(uv|ii) for all u&v

        for ( int t_c = first_index_[h][1]; t_c < last_index_[h][1]; t_c++){

            int tt_den = Active_GindMap_c_[Full_ij_index(t_c,t_c,nmo_)];

            for ( int a_c = first_index_[h][2]; a_c < last_index_[h][2]; a_c++){

                int aa = a_c - first_index_[h][2];

                double Hval = 0.0e0;

                for ( int u_c = first_index_[h_u][1]; u_c < last_index_[h_u][1]; u_c++){

                    // d2(tt|uv)*tei(uv|aa) ... (u > v)

                    for ( int v_c = first_index_[h_u][1]; v_c < u_c; v_c++){

                        int uv_den   = Active_GindMap_c_[Full_ij_index(u_c,v_c,nmo_)];
                        int ttuv_den = Lt_ij_index(tt_den,uv_den);
                        int aauv_int = Full_ij_index(uv_den,aa,ngem_ts);

                        Hval += d2[ttuv_den] * scr_4index_1[aauv_int];

                    }

                    // 0.5 * d2(tt|uu)*tei(uu|aa)      

                    int uu_den   = Active_GindMap_c_[Full_ij_index(u_c,u_c,nmo_)];
                    int ttuu_den = Lt_ij_index(tt_den,uu_den);
                    int aauu_int = Full_ij_index(uu_den,aa,ngem_ts);

                    Hval += 0.5e0 * d2[ttuu_den] * scr_4index_1[aauu_int];

                }

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(t_c,a_c)];

                // 2*d2(tt|uv)*tei(uv|aa) ... (for all u/v)

                diagonal_orbital_hessian_[grad_ind] += 4.0e0 * Hval;

            }

        }

    }

    // exchange terms

    for ( int a_c = first_index_[h][2]; a_c < last_index_[h][2]; a_c++ ){

        for ( int h_u = 0; h_u < nirrep_; h_u++) {

            int amo_offset = first_index_[h_u][1];
            int nmo_h_u    = amopi_[h_u];

            if ( nmo_h_u == 0 ) continue;

            for ( int u_c = first_index_[h_u][1]; u_c < last_index_[h_u][1]; u_c++){

                int ua_df  = df_ij_index(u_c,a_c);
                int ua_scr = ( u_c - amo_offset ) * nQ_;

                C_DCOPY(nQ_,&tei[ua_df],Qstride_,&scr_3index_1[ua_scr],1);

            }

            // calculate all (ua|va) for this a and u/v in h_u

            F_DGEMM('t','n',nmo_h_u,nmo_h_u,nQ_,1.0e0,&scr_3index_1[0],nQ_,
                   &scr_3index_1[0],nQ_,0.0e0,&scr_4index_1[0],nmo_h_u);

            // update Hessian elements

            for ( int t_c = first_index_[h][1]; t_c < last_index_[h][1]; t_c++){

                double Hval = 0.0e0;

                for ( int u_c = first_index_[h_u][1]; u_c < last_index_[h_u][1]; u_c++){

                    // d2(tu|tv)*tei(au|av) ... ( u > v )

                    for ( int v_c = first_index_[h_u][1]; v_c < u_c; v_c++){

                        int ua_scr   = ( u_c - amo_offset );
                        int va_scr   = ( v_c - amo_offset );
                        int uava_int = Full_ij_index(ua_scr,va_scr,nmo_h_u);

                        int ut_den   = Active_GindMap_c_[Full_ij_index(u_c,t_c,nmo_)];
                        int vt_den   = Active_GindMap_c_[Full_ij_index(v_c,t_c,nmo_)];
                        int utvt_den = Lt_ij_index(ut_den,vt_den);

                        Hval += scr_4index_1[uava_int] * d2[utvt_den];

                    }

                    // 0.5 * d2(tu|tu)*tei(au|au)

                    int ua_scr   = ( u_c - amo_offset );
                    int uaua_int = Full_ij_index(ua_scr,ua_scr,nmo_h_u);

                    int ut_den   = Active_GindMap_c_[Full_ij_index(u_c,t_c,nmo_)];
                    int utut_den = Lt_ij_index(ut_den,ut_den);

                    Hval += 0.5e0 * scr_4index_1[uaua_int] * d2[utut_den];

                }

                int grad_ind = Grad_IndMap_c_[Lt_ij_index(t_c,a_c)];

                // 4.0 * d2(tu|tv)*tei(au|av) ... (for all u/v)

                diagonal_orbital_hessian_[grad_ind] += 8.0e0 * Hval;

            }

        }

    }

}

void OrbitalOptimizer::Alloc_Hessian(){

    int Nrot = 0;

    for ( int p_class = 0; p_class < 3; p_class++ ){

        int q_class_max = p_class;
        if ( active_active_rotations_ && p_class == 1 ) q_class_max = p_class + 1;

        for ( int q_class = 0; q_class < q_class_max; q_class++ ){

            for ( int h = 0; h < nirrep_; h++ ){

                int nmo_p = last_index_[h][p_class] - first_index_[h][p_class];

                int nmo_q = last_index_[h][q_class] - first_index_[h][q_class];

                if ( p_class == q_class ){

                   Nrot += nmo_p * ( nmo_p - 1 ) / 2;

                }

                else {

                   Nrot += nmo_p * nmo_q;

                }

            } // end h loop

        } // emd q_class loop       

    } // end p_class loop

    if ( Nrot_ != Nrot ) {

        throw PsiException("Nrot dimension mismatch when allocating Hessian",__FILE__,__LINE__);

    }

    diagonal_orbital_hessian_ = (double *)malloc(Nrot_*sizeof(double));
    memset((void*)diagonal_orbital_hessian_,'\0',Nrot_*sizeof(double));

}

void OrbitalOptimizer::Dealloc_Hessian(){

    free(diagonal_orbital_hessian_);

}




}
