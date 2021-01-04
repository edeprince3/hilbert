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

void OrbitalOptimizer::TransformIntegrals( double * tei, double * oei, double * MOcoeff ){

    TransformOei( oei );

    TransformTei( tei );

    TransformMOcoeff( MOcoeff );

}

void OrbitalOptimizer::TransformMOcoeff(double * MOcoeff){

    for ( int h = 0; h < nirrep_; h++){

        int nmo_h = nmopi_[h];
 
        if ( nmo_h == 0 ) continue;

        C_DCOPY(nmo_h*nmo_h,&MOcoeff[U_offset_[h]],1,&Tei_scr1_[0],1);

        F_DGEMM('n','n',nmo_h,nmo_h,nmo_h,1.0e0,&Tei_scr1_[0],nmo_h,\
                                                &U_[U_offset_[h]],nmo_h, \
                                          0.0e0,&Tei_scr2_[0],nmo_h);

        F_DGEMM('t','n',nmo_h,nmo_h,nmo_h,1.0e0,&U_[U_offset_[h]],nmo_h, \
                                                &Tei_scr2_[0],nmo_h, \
                                          0.0e0,&Tei_scr1_[0],nmo_h);

        C_DCOPY(nmo_h*nmo_h,&Tei_scr1_[0],1,&MOcoeff[U_offset_[h]],1);

    }


}

void OrbitalOptimizer::TransformTei(double * tei){

    for ( int Q = 0; Q < nQ_; Q++){

        int thread_id = 0;

        TransformTei_Q(&tei[Q*Qstride_], \
                       &Tei_scr1_[thread_id*Tei_scr_dim_], \
                       &Tei_scr2_[thread_id*Tei_scr_dim_]);

    }

}

void OrbitalOptimizer::Print_IntQ(double * tei){

    for ( int h_p = 0; h_p < nirrep_; h_p++){

        for ( int h_q = 0; h_q < h_p + 1; h_q++){

            for ( int p_class = 0; p_class < 3; p_class++){

                for ( int q_class = 0; q_class < 3; q_class++){

                    for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                        int q_c_max = last_index_[h_q][q_class];
                        if ( q_c_max > last_index_[h_p][p_class] ) q_c_max = p_c + 1;

                        for ( int q_c = first_index_[h_q][q_class]; q_c < q_c_max; q_c++){

                            printf("%i %i    %i %i    %15.7e\n",h_p,h_q,p_c,q_c,tei[df_ij_index(p_c,q_c)]);

                        }

                    }

                }

            }

        }

    }

}

void OrbitalOptimizer::TransformTei_UnpackQ(double * teiQ, \
                                            double * tei_scr1,
                                            int h_p, int h_q){

    int nmo_p = nmopi_[h_p];
    int nmo_q = nmopi_[h_q];

    if ( h_p == h_q ){

        for ( int q_class = 0; q_class < 3; q_class++){

            for ( int q_c = first_index_[h_q][q_class]; q_c < last_index_[h_q][q_class]; q_c++){

                int q_i       = OindMap_c2i_[q_c];

                for ( int p_class = 0; p_class < 3; p_class++){

                    for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                        int p_i       = OindMap_c2i_[p_c];
                        int tei_ind   = df_ij_index(p_c,q_c);

                        tei_scr1[Full_ij_index(p_i,q_i,nmo_p)] = teiQ[tei_ind];

                    }

                }

            }

        }

    }

    else{

        for ( int q_class = 0; q_class < 3; q_class++){

            for ( int q_c = first_index_[h_q][q_class]; q_c < last_index_[h_q][q_class]; q_c++){

                int q_i       = OindMap_c2i_[q_c];

                for ( int p_class = 0; p_class < 3; p_class++){

                    for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                        int p_i       = OindMap_c2i_[p_c];
                        int tei_ind   = df_ij_index(p_c,q_c);

                        tei_scr1[Full_ij_index(p_i,q_i,nmo_p)] = teiQ[tei_ind];

                    }

                }

            }

        }

    }


}

void OrbitalOptimizer::TransformTei_PackQ(double * teiQ, \
                                            double * tei_scr1,
                                            int h_p, int h_q){

    int nmo_p = nmopi_[h_p];
    int nmo_q = nmopi_[h_q];

    if ( h_p == h_q ){

        for ( int q_class = 0; q_class < 3; q_class++){

            for ( int q_c = first_index_[h_q][q_class]; q_c < last_index_[h_q][q_class]; q_c++){

                int q_i       = OindMap_c2i_[q_c];

                for ( int p_class = 0; p_class < 3; p_class++){

                    for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                        int p_i       = OindMap_c2i_[p_c];
                        int tei_ind   = df_ij_index(p_c,q_c);

                        teiQ[tei_ind] = tei_scr1[Full_ij_index(p_i,q_i,nmo_p)];

                    }

                }

            }

        }

    }

    else{


        for ( int q_class = 0; q_class < 3; q_class++){

            for ( int q_c = first_index_[h_q][q_class]; q_c < last_index_[h_q][q_class]; q_c++){

                int q_i       = OindMap_c2i_[q_c];

                for ( int p_class = 0; p_class < 3; p_class++){

                    for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                        int p_i       = OindMap_c2i_[p_c];
                        int tei_ind   = df_ij_index(p_c,q_c);

                        teiQ[tei_ind] = tei_scr1[Full_ij_index(p_i,q_i,nmo_p)];

                    }

                }

            }

        }

    }

}

void OrbitalOptimizer::TransformTei_Q(double * teiQ, double * tei_scr1, double * tei_scr2){

    for ( int h_p = 0; h_p < nirrep_; h_p++ ){

        for ( int h_q = 0; h_q < h_p + 1; h_q++ ){
            
            int nmo_p = nmopi_[h_p];
            int nmo_q = nmopi_[h_q];

            if ( nmo_p == 0 || nmo_q == 0 ) continue;

            TransformTei_UnpackQ(&teiQ[0], &tei_scr1[0], h_p, h_q);

            F_DGEMM('n','n',nmo_p,nmo_q,nmo_q,1.0e0,&tei_scr1[0],nmo_p,        \
                                                    &U_[U_offset_[h_q]],nmo_q, \
                                              0.0e0,&tei_scr2[0],nmo_p);

            F_DGEMM('t','n',nmo_p,nmo_q,nmo_p,1.0e0,&U_[U_offset_[h_p]],nmo_p, \
                                                    &tei_scr2[0],nmo_p,\
                                              0.0e0,&tei_scr1[0],nmo_p);

            TransformTei_PackQ(&teiQ[0], &tei_scr1[0], h_p, h_q);

        }

    }
}

void OrbitalOptimizer::TransformOei(double * oei){

   for ( int h = 0; h < nirrep_; h++ ){

       int nmo_p = nmopi_[h];

       if ( nmo_p == 0 ) continue;

       TransformOei_UnpackSymBlock( oei , h);

       F_DGEMM('n','n',nmo_p,nmo_p,nmo_p,1.0e0,&Tei_scr1_[0],nmo_p,     \
                                               &U_[U_offset_[h]],nmo_p, \
                                         0.0e0,&Tei_scr2_[0],nmo_p);

       F_DGEMM('t','n',nmo_p,nmo_p,nmo_p,1.0e0,&U_[U_offset_[h]],nmo_p, \
                                                   &Tei_scr2_[0],nmo_p, \
                                         0.0e0,&Tei_scr1_[0],nmo_p);

       TransformOei_PackSymBlock( oei , h);

   }

}

void OrbitalOptimizer::TransformOei_PackSymBlock(double * oei, int h ){

   int nmo_p = nmopi_[h];

   for ( int p_class = 0; p_class < 3; p_class++ ){

       for ( int p_c = first_index_[h][p_class]; p_c < last_index_[h][p_class]; p_c++ ){

           int p_i = OindMap_c2i_[p_c];

           for ( int q_c = p_c; q_c < last_index_[h][p_class]; q_c++ ){

               int q_i       = OindMap_c2i_[q_c];

               int oei_index = Full_GindMap_c_[Full_ij_index(p_c,q_c,nmo_)];

               int scr_index = Full_ij_index(p_i,q_i,nmo_p);

               oei[oei_index] = Tei_scr1_[Full_ij_index(q_i,p_i,nmo_p)];

           }

           for ( int q_class = p_class + 1; q_class < 3; q_class++ ){

               for ( int q_c = first_index_[h][q_class]; q_c < last_index_[h][q_class]; q_c++ ){

                   int q_i       = OindMap_c2i_[q_c];

                   int oei_index = Full_GindMap_c_[Full_ij_index(p_c,q_c,nmo_)];

                   oei[oei_index] = Tei_scr1_[Full_ij_index(q_i,p_i,nmo_p)];

               }

           }

       }

   }

}

void OrbitalOptimizer::TransformOei_UnpackSymBlock(double * oei, int h ){

   int nmo_p = nmopi_[h];

   for ( int p_class = 0; p_class < 3; p_class++ ){

       for ( int p_c = first_index_[h][p_class]; p_c < last_index_[h][p_class]; p_c++ ){

           int p_i = OindMap_c2i_[p_c];

           for ( int q_c = p_c; q_c < last_index_[h][p_class]; q_c++ ){

               int q_i       = OindMap_c2i_[q_c];

               int oei_index = Full_GindMap_c_[Full_ij_index(p_c,q_c,nmo_)];

               int scr_index = Full_ij_index(p_i,q_i,nmo_p);

               Tei_scr1_[Full_ij_index(p_i,q_i,nmo_p)] = oei[oei_index];
               Tei_scr1_[Full_ij_index(q_i,p_i,nmo_p)] = oei[oei_index];

           }

           for ( int q_class = p_class + 1; q_class < 3; q_class++ ){

               for ( int q_c = first_index_[h][q_class]; q_c < last_index_[h][q_class]; q_c++ ){

                   int q_i       = OindMap_c2i_[q_c];

                   int oei_index = Full_GindMap_c_[Full_ij_index(p_c,q_c,nmo_)];

                   Tei_scr1_[Full_ij_index(p_i,q_i,nmo_p)] = oei[oei_index];
                   Tei_scr1_[Full_ij_index(q_i,p_i,nmo_p)] = oei[oei_index];

               }

           }

       }

   }

}




}
