#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libqt/qt.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libtrans/mospace.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>

#include <misc/blas.h>
#include <misc/omp.h>

#include "../orbital_optimizer.h"

using namespace psi;
using namespace fnocc;

namespace hilbert{


double OrbitalOptimizer::ComputeE(double * d2, double * d1, double * tei, double * oei){

    // function o calculate the electronic energy

    E1_c_  = ComputeE1_c(oei);

    E1_a_  = ComputeE1_a(d1, oei);

    E2_cc_ = ComputeE2_cc(tei);

    E2_ca_ = ComputeE2_ca(d1, tei);

    E2_aa_ = ComputeE2_aa(d2, tei);

    E_elec_ = E1_c_ + E1_a_ + E2_cc_ + E2_ca_ + E2_aa_;

    return E_elec_;

}

double OrbitalOptimizer::ComputeE_2index(double * d1, double * oei){

    // function to calculate energy using Fi_, Fa_, Q_, and Z_, oei, and d1
    // E = 0.5 * sum_{pq} d1[pq]*oei[pq]
    //         + sum{p}   F[pp]    
    // where F is the generalized Fock matrix
    //
    // E = sum_{i} oei[ii] + sum_{i} Fi_[ii] + sum_{i} Fa_[ii]
    //     0.5 * sum_{uv} d1[uv] * oei[uv] + 
    //     0.5 * sum_{u} Z_[uu] + 0.5 * sum_{u} Q_[uu]

    double E = 0.0e0;

    for ( int h = 0; h < nirrep_; h++ ){

        for ( int i_c = first_index_[h][0]; i_c < last_index_[h][0]; i_c++ ){

            E += oei[ Full_GindMap_c_[Full_ij_index(i_c,i_c,nmo_)] ]
               + Fi_[ Lt_ij_oe_index(i_c,i_c,h) ]
               + Fa_[ Lt_ij_oe_index(i_c,i_c,h) ];

        }

        for ( int u_c = first_index_[h][1]; u_c < last_index_[h][1]; u_c++ ){

            E += 0.5e0 * (
               oei[ Full_GindMap_c_[Full_ij_index(u_c,u_c,nmo_)] ]
             * d1[ Active_GindMap_c_[Full_ij_index(u_c,u_c,nmo_)] ]
             + Z_[ Q_index(u_c,u_c,h) ]
             + Q_[ Q_index(u_c,u_c,h) ]
             );

             for ( int v_c = first_index_[h][1]; v_c < u_c; v_c++ ){

                 E += oei[ Full_GindMap_c_[Full_ij_index(u_c,v_c,nmo_)] ]
                    *  d1[ Active_GindMap_c_[Full_ij_index(u_c,v_c,nmo_)] ];

             }

        }

    }

    E_elec_ = E;

    return E;

}

double OrbitalOptimizer::ComputeE2_aa(double * d2, double * tei){

    // function calculate 2-e contribution to energy (active-active)

    double E2_aa=0.0e0;

    int max_gempi = Max_value(act_gempi_,nirrep_);

    double * scr_3index_1 = (double *)malloc(max_gempi*nQ_*sizeof(double));
    memset((void*)scr_3index_1,'\0',max_gempi*nQ_*sizeof(double));

    double * scr_4index_1 = (double *)malloc(max_gempi*max_gempi*sizeof(double));
    memset((void*)scr_4index_1,'\0',max_gempi*max_gempi*sizeof(double));

    int offset = 0;

    for ( int h_uv = 0; h_uv < nirrep_; h_uv++){

        int ngem_uv = act_gempi_[h_uv];

        if ( ngem_uv == 0 ) continue;

        // gather all (:,uv) active-active

        for ( int h_u = 0; h_u < nirrep_; h_u++){

           int h_v = GemSym(h_u,h_uv);

           if ( h_v > h_u ) continue;

           //$pragma omp parallel for schedule (dynamic)
           for ( int u_c = first_index_[h_u][1]; u_c < last_index_[h_u][1]; u_c++){

               int max_v_c = last_index_[h_v][1];
               if ( h_v == h_u ) max_v_c = u_c + 1;

               for ( int v_c = first_index_[h_v][1]; v_c < max_v_c; v_c++ ){

                   int uv_df  = df_ij_index(u_c,v_c);
                   int uv_den = Active_GindMap_c_[Full_ij_index(u_c,v_c,nmo_)];

                   C_DCOPY(nQ_,&tei[uv_df],Qstride_,&scr_3index_1[uv_den*nQ_],1);

               } // end v_c loop

           } // end u_c loop

        } // end h_u loop

        // calculate all (uv|wx) integrals

        F_DGEMM('t','n',ngem_uv,ngem_uv,nQ_,1.0e0,&scr_3index_1[0],nQ_,&scr_3index_1[0],nQ_,0.0e0,&scr_4index_1[0],ngem_uv);

        for ( int h_u = 0; h_u < nirrep_; h_u++){

            int h_v = GemSym(h_uv,h_u);

            for ( int h_w = 0; h_w < nirrep_; h_w++){

                int h_x = GemSym(h_uv,h_w);

                for ( int u_c = first_index_[h_u][1]; u_c < last_index_[h_u][1]; u_c++){

                    for ( int v_c = first_index_[h_v][1]; v_c < last_index_[h_v][1]; v_c++){

                        int uv_den = Active_GindMap_c_[Full_ij_index(u_c,v_c,nmo_)];

                        for ( int w_c = first_index_[h_w][1]; w_c < last_index_[h_w][1]; w_c++){

                            for ( int x_c = first_index_[h_x][1]; x_c < last_index_[h_x][1]; x_c++){

                                int wx_den   = Active_GindMap_c_[Full_ij_index(w_c,x_c,nmo_)];
                                int uvwx_int = Full_ij_index(uv_den,wx_den,ngem_uv);
                                int uvwx_den = Lt_ij_index(uv_den,wx_den) + offset;

                                E2_aa += d2[uvwx_den] * scr_4index_1[uvwx_int];

                            } // end u_c loop

                        } // end v_c loop


                    } // end w_c loop


                } // end x_c loop

                

            } // end h_w loop

        } // end h_u loop

        offset += ngem_uv * ( ngem_uv + 1 ) / 2;

    } // end h_uv loop

    free(scr_3index_1);
    free(scr_4index_1);

    return 5.0e-1 * E2_aa;

}

double OrbitalOptimizer::ComputeE2_ca(double * d1, double * tei){

    // function to calculate 2-e contribution to energy (core-active)

    double E2_ca = 0.0e0;

    int max_gempi = Max_value(act_gempi_,nirrep_);

    int * scr_offset = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)scr_offset,'\0',nirrep_*sizeof(int));

    double * scr_3index_1 = (double *)malloc(nrstc_*nQ_*sizeof(double));
    memset((void*)scr_3index_1,'\0',nrstc_*nQ_*sizeof(double));

    double * scr_3index_2 = (double *)malloc(nrstc_*nQ_*sizeof(double));
    memset((void*)scr_3index_2,'\0',nrstc_*nQ_*sizeof(double));

    double * scr_4index_1 = (double *)malloc(max_gempi*sizeof(double));
    memset((void*)scr_4index_1,'\0',max_gempi*sizeof(double));

    double * scr_4index_2 = (double *)malloc(max_gempi*sizeof(double));
    memset((void*)scr_4index_2,'\0',max_gempi*sizeof(double));

    double * vec_uv = (double *)malloc(nQ_*sizeof(double));
    memset((void*)vec_uv,'\0',nQ_*sizeof(double));

    double * vec_uu = (double *)malloc(nQ_*sizeof(double));
    memset((void*)vec_uu,'\0',nQ_*sizeof(double));

    for (int h_i = 1; h_i < nirrep_; h_i++) scr_offset[h_i] = scr_offset[h_i - 1 ] + rstcpi_[h_i - 1];

    for ( int h_i = 0; h_i < nirrep_; h_i++){

        //* gather all (:,ii) integrals core-core

        //$pragma omp parallel for schedule (static)
        for ( int i_c = first_index_[h_i][0]; i_c < last_index_[h_i][0]; i_c++){

            int ii_df = df_ij_index(i_c,i_c);
            int i_i   = OindMap_c2i_[i_c] + scr_offset[h_i];

            C_DCOPY(nQ_,&tei[ii_df],Qstride_,&scr_3index_1[i_i*nQ_],1);

        } // end i_c loop

    } // end h_i loop

    for ( int h_i = 0; h_i < nirrep_; h_i++){

        for ( int u_c = first_index_[h_i][1]; u_c < last_index_[h_i][1]; u_c++){

            //* gather all (:,uj) integrals active-core

            for ( int h_k = 0; h_k < nirrep_; h_k++){

                //$pragma omp parallel for schedule (static)
                for ( int k_c = first_index_[h_k][0]; k_c < last_index_[h_k][0]; k_c++){

                    int uk_df = df_ij_index(u_c,k_c);
                    int k_i   = OindMap_c2i_[k_c] + scr_offset[h_k];

                    C_DCOPY(nQ_,&tei[uk_df],Qstride_,&scr_3index_2[k_i*nQ_],1);

                } // end k_c loop

            } // end h_k loop

            for ( int v_c = first_index_[h_i][1]; v_c < u_c; v_c++){

                int uv_den = Active_GindMap_c_[Full_ij_index(u_c,v_c,nmo_)];
                int uv_df  = df_ij_index(u_c,v_c);

                C_DCOPY(nQ_,&tei[uv_df],Qstride_,&vec_uv[0],1);

                for ( int h_k = 0; h_k < nirrep_; h_k++){

                    for ( int k_c = first_index_[h_k][0]; k_c < last_index_[h_k][0]; k_c++){

                         int k_i   = OindMap_c2i_[k_c] + scr_offset[h_k];
                         int vk_df = df_ij_index(v_c,k_c);

                         scr_4index_1[uv_den] += 2.0e0 * C_DDOT(nQ_,&vec_uv[0],1,&scr_3index_1[k_i*nQ_],1);
                         scr_4index_2[uv_den] += C_DDOT(nQ_,&tei[vk_df],Qstride_,&scr_3index_2[k_i*nQ_],1);

                    } // end k_c loop

                } // end h_k looop                

            } // end v_c loop

            int uu_den = Active_GindMap_c_[Full_ij_index(u_c,u_c,nmo_)];
            int uu_df  = df_ij_index(u_c,u_c);

            C_DCOPY(nQ_,&tei[uu_df],Qstride_,&vec_uu[0],1);

            for ( int h_k = 0; h_k < nirrep_; h_k++){

                for ( int k_c = first_index_[h_k][0]; k_c < last_index_[h_k][0]; k_c++){

                     int k_i   = OindMap_c2i_[k_c] + scr_offset[h_k];

                     scr_4index_1[uu_den] += C_DDOT(nQ_,&vec_uu[0],1,&scr_3index_1[k_i*nQ_],1);
                     scr_4index_2[uu_den] += 5.0e-1 * C_DDOT(nQ_,&scr_3index_2[k_i*nQ_],1,&scr_3index_2[k_i*nQ_],1);

                 } // end k_c loop

            } // end h_k loop

        } // end u_c for

    } // end h_i for

    for ( int uv = 0; uv < max_gempi; uv++) E2_ca += ( scr_4index_1[uv] - scr_4index_2[uv] ) * d1[uv];

    E2_ca = 2.0e0 * E2_ca;

    free(scr_3index_1);
    free(scr_3index_2);
    free(scr_4index_1);
    free(scr_4index_2);
    free(vec_uv);
    free(vec_uu);
    free(scr_offset);

    return E2_ca;

}

double OrbitalOptimizer::ComputeE2_cc(double * tei){

    // function to calculate 2-e contribution to energy (core-core)

    double E2_cc = 0.0e0;

    int max_docpi = Max_value(rstcpi_,nirrep_);

    double * scr_3index_1 = (double *)malloc(max_docpi*nQ_*sizeof(double));
    memset((void*)scr_3index_1,'\0',max_docpi*nQ_*sizeof(double));

    double * vec_tmp = (double *)malloc(nQ_*sizeof(double));
    memset((void*)vec_tmp,'\0',nQ_*sizeof(double));

    double E_J = 0.0e0;
    double E_K = 0.0e0;

    for ( int h_i = 0; h_i < nirrep_; h_i++){

        //$pragma omp parallel for schedule (static)
        for ( int i_c = first_index_[h_i][0]; i_c < last_index_[h_i][0]; i_c++){

            int ii_df = df_ij_index(i_c,i_c);
            int i_i   = OindMap_c2i_[i_c];

            C_DCOPY(nQ_,&tei[ii_df],Qstride_,&scr_3index_1[i_i*nQ_],1);

        }

        for ( int h_j = 0; h_j < h_i; h_j++){

            for ( int i_c = first_index_[h_i][0]; i_c < last_index_[h_i][0]; i_c++){

                int i_i  = OindMap_c2i_[i_c];

                for ( int j_c = first_index_[h_j][0]; j_c < last_index_[h_j][0]; j_c++){

                    int ij_df = df_ij_index(i_c,j_c);
                    int jj_df = df_ij_index(j_c,j_c);

                    E_J += C_DDOT(nQ_,&scr_3index_1[i_i*nQ_],1,&tei[jj_df],Qstride_);
                    E_K += C_DDOT(nQ_,&tei[ij_df],Qstride_,&tei[ij_df],Qstride_);

                }

            }

        }

        for ( int i_c = first_index_[h_i][0]; i_c < last_index_[h_i][0]; i_c++){

            int i_i  = OindMap_c2i_[i_c];

            double tmp = 5.0e-1 * C_DDOT(nQ_,&scr_3index_1[i_i*nQ_],1,&scr_3index_1[i_i*nQ_],1);
            E_J += tmp;
            E_K += tmp;

            for ( int j_c = first_index_[h_i][0]; j_c < i_c; j_c++){

                int j_i   = OindMap_c2i_[j_c];
                int ij_df = df_ij_index(i_c,j_c);

                E_J += C_DDOT(nQ_,&scr_3index_1[i_i*nQ_],1,&scr_3index_1[j_i*nQ_],1);
                E_K += C_DDOT(nQ_,&tei[ij_df],Qstride_,&tei[ij_df],Qstride_);

            }

        }

    }

    free(scr_3index_1);

    E2_cc = 2.0e0 * ( 2.0e0 * E_J - E_K );

    return E2_cc;

}

double OrbitalOptimizer::ComputeE1_a(double * d1, double * oei){

    // function to calculate 1-e contribution to energy (active)

    double E1_a =0.0e0;

    for ( int h_p = 0; h_p < nirrep_; h_p++){

        for ( int p_c = first_index_[h_p][1]; p_c < last_index_[h_p][1]; p_c++){

            for ( int q_c = first_index_[h_p][1]; q_c < p_c; q_c++){

                int pq_d1  = Active_GindMap_c_[Full_ij_index(p_c,q_c,nmo_)];
                int pq_oei = Full_GindMap_c_[Full_ij_index(p_c,q_c,nmo_)];

                E1_a += 2.0e0 * d1[pq_d1] * oei[pq_oei];

            }

            int pp_d1  = Active_GindMap_c_[Full_ij_index(p_c,p_c,nmo_)];
            int pp_oei = Full_GindMap_c_[Full_ij_index(p_c,p_c,nmo_)];

            E1_a += d1[pp_d1] * oei[pp_oei];

        }

    }

     return E1_a;

}

double OrbitalOptimizer::ComputeE1_c(double * oei){

    // function to calculate 1-e contribution to energy (core)

    double E1_c =0.0e0;

    for ( int h_p = 0; h_p < nirrep_; h_p++){

        for ( int p_c = first_index_[h_p][0]; p_c < last_index_[h_p][0]; p_c++){

            int pp_oei = Full_GindMap_c_[Full_ij_index(p_c,p_c,nmo_)];

            E1_c += 2.0e0 * oei[pp_oei];

        }

    }

    return E1_c;

}



}
