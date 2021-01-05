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

void OrbitalOptimizer::Sort(double * d2, double * d1, double * oei, int direction){

    // function to sort d2, d1, and oei from energy to class order (direction 1) or class to energy (direction -1)
    // elements of d2 are scaled as well

    int d2_max_ngem  = Max_value(act_gempi_,nirrep_);
    int oei_max_ngem = full_gempi_[0];

    int tmp_dim = d2_max_ngem*(d2_max_ngem+1)/2;
    if ( oei_max_ngem > tmp_dim ) tmp_dim = oei_max_ngem;

    double * symblock_tmp = (double *)malloc(tmp_dim*sizeof(double));
    memset((void*)symblock_tmp,'\0',tmp_dim*sizeof(double));

    // sort 2-e density

    int irrep_offset = 0;

    for ( int h_pq = 0; h_pq < nirrep_; h_pq++){

        int block_dim = act_gempi_[h_pq]*(act_gempi_[h_pq]+1)/2;

        // sort and scale density elements and save into vector

        for ( int h_p = 0; h_p < nirrep_; h_p++){

            int h_q = GemSym(h_p,h_pq);

            if ( h_q > h_p ) continue;

            for ( int h_r = 0; h_r < nirrep_; h_r++){

                int h_s = GemSym(h_r,h_pq);

                if ( h_s > h_r ) continue;

                for ( int p_c = first_index_[h_p][1]; p_c < last_index_[h_p][1]; p_c++){

                    int max_q_c = last_index_[h_q][1];
                    if ( h_p == h_q ) max_q_c = p_c + 1;

                    int p_e = OindMap_c2e_[p_c];

                    for ( int q_c = first_index_[h_q][1]; q_c < max_q_c; q_c++ ) {

                        double pq_fac = 5.0e-1;
                        if ( p_c == q_c ) pq_fac = 1.0e0;

                        int q_e = OindMap_c2e_[q_c];

                        int pq_c = Active_GindMap_c_[Full_ij_index(p_c,q_c,nmo_)];
                        int pq_e = Active_GindMap_e_[Full_ij_index(p_e,q_e,nmo_)];

                        for ( int r_c = first_index_[h_r][1]; r_c < last_index_[h_r][1]; r_c++){

                            int max_s_c = last_index_[h_s][1];
                            if ( h_r == h_s ) max_s_c = r_c + 1;

                            int r_e = OindMap_c2e_[r_c];

                            for ( int s_c = first_index_[h_s][1]; s_c < max_s_c; s_c++ ) {

                                double rs_fac = 5.0e-1;
                                if ( r_c == s_c ) rs_fac = 1.0e0;

                                int s_e = OindMap_c2e_[s_c];

                                int rs_c = Active_GindMap_c_[Full_ij_index(r_c,s_c,nmo_)];
                                int rs_e = Active_GindMap_e_[Full_ij_index(r_e,s_e,nmo_)];

                                double pqrs_fac = 5.0e-1*pq_fac*rs_fac;
                                if ( pq_c == rs_c ) pqrs_fac = pq_fac*rs_fac;

                                if ( direction == 1  ) symblock_tmp[Lt_ij_index(pq_c,rs_c)] = pqrs_fac * d2[irrep_offset + Lt_ij_index(pq_e,rs_e)];
                                if ( direction == -1 ) symblock_tmp[Lt_ij_index(pq_e,rs_e)] = d2[irrep_offset + Lt_ij_index(pq_c,rs_c)]/pqrs_fac;

                            }

                        }

                    }

                }

            }

        }

        // copy reordered block into old one

        if ( direction == 1 ) for (int i = 0; i < block_dim; i++) d2[irrep_offset + i] = 2.0e0 * symblock_tmp[i];

        if ( direction == -1 ) for (int i = 0; i < block_dim; i++) d2[irrep_offset + i] = 5.0e-1 * symblock_tmp[i];

        irrep_offset += block_dim;
    }

    // sort 1-e density

    irrep_offset = 0;

    for ( int h_p = 0; h_p < nirrep_; h_p++ ){

        for ( int p_c = first_index_[h_p][1]; p_c < last_index_[h_p][1]; p_c++){

            int p_i = OindMap_c2i_[p_c] - rstcpi_[h_p];

            for ( int q_c = first_index_[h_p][1]; q_c < p_c + 1; q_c++){

                int q_i = OindMap_c2i_[q_c] - rstcpi_[h_p];

                int pq_i = Lt_ij_index(p_i,q_i)+irrep_offset;

                int pq_c = Active_GindMap_c_[Full_ij_index(p_c,q_c,nmo_)];

                if ( direction == 1  ) symblock_tmp[pq_c] = d1[pq_i];
                if ( direction == -1 ) symblock_tmp[pq_i] = d1[pq_c];

            }

        }

        irrep_offset += amopi_[h_p] * ( amopi_[h_p] + 1 ) / 2;

    }


    for ( int i = 0; i < act_gempi_[0]; i++) d1[i] = symblock_tmp[i];

    // sort 1-e integrals

    irrep_offset = 0;

    for ( int h_p = 0; h_p < nirrep_; h_p++){

        for (int p_class = 0; p_class < 3; p_class++){

            for ( int q_class = 0; q_class <= p_class; q_class++){

                for (int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                   int p_i = OindMap_c2i_[p_c];

                   int q_c_max = last_index_[h_p][q_class];
                   if ( p_class == q_class ) q_c_max = p_c + 1;

                   for ( int q_c = first_index_[h_p][q_class]; q_c < q_c_max; q_c++){

                       int q_i = OindMap_c2i_[q_c];

                       int pq_c = Full_GindMap_c_[Full_ij_index(p_c,q_c,nmo_)];
                       int pq_i = Lt_ij_index(p_i,q_i)+irrep_offset;

                       if ( direction == 1  ) symblock_tmp[pq_c] = oei[pq_i];
                       if ( direction == -1 ) symblock_tmp[pq_i] = oei[pq_c];

                   }

                }

            }

        }

        irrep_offset += nmopi_[h_p] * ( nmopi_[h_p] + 1 ) / 2;

    }

    for ( int i = 0; i < full_gempi_[0]; i++) oei[i] = symblock_tmp[i];

    free(symblock_tmp);

}


void OrbitalOptimizer::Sort_MOcoeff(double * trans, int direction){

    // sort the MO coefficient matrix

    double * symblock_tmp;

    if ( direction == -1 ){

        symblock_tmp = (double *)malloc(nmo_*nmo_*sizeof(double));
        memset((void*)symblock_tmp,'\0',nmo_*nmo_*sizeof(double));

    }

    int irrep_offset = 0;

    for ( int h = 0; h < nirrep_; h++ ){

        int nmo = nmopi_[h]; 

        for (int p_class = 0; p_class < 3; p_class++){

            for ( int p_c = p_c = first_index_[h][p_class]; p_c < last_index_[h][p_class]; p_c++){

                int p_i = OindMap_c2i_[p_c];
                int p_e = OindMap_c2e_[p_c];

                for (int q_class = 0; q_class < 3; q_class++){

                    for ( int q_c = first_index_[h][q_class]; q_c < last_index_[h][q_class]; q_c++){
 
                        int q_i  = OindMap_c2i_[q_c];
                        int q_e  = OindMap_c2e_[q_c];
                        int pq_i = Full_ij_index(q_i,p_i,nmo) + irrep_offset;
                        int pq_e = Full_ij_index(q_e,p_e,nmo_);

                        if ( direction == 1 ) {

                            if ( p_i == q_i ) {

                                trans[pq_i] = 1.0e0;

                            }

                            else {

                                trans[pq_i] = 0.0e0;
       
                            }

                        }

                        else {

                            symblock_tmp[pq_e] = trans[pq_i];

                        }
  
                    }

                }

            } 

        }

        irrep_offset += nmo * nmo;

    }

    if ( direction == -1 ) {

        C_DCOPY(nmo_*nmo_,&symblock_tmp[0],1,&trans[0],1);
        free(symblock_tmp);
 
    }

}



}
