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

#include "../orbital_optimizer.h"

using namespace psi;
using namespace fnocc;

namespace hilbert{

void OrbitalOptimizer::Determine_OindMap(){

    Qstride_ = nmo_ * ( nmo_ + 1 ) / 2;

    // Symmetry direct product table

    SymProd_[0*8+1] = SymProd_[1*8+0] = 1;
    SymProd_[0*8+2] = SymProd_[2*8+0] = 2;
    SymProd_[0*8+3] = SymProd_[3*8+0] = 3;
    SymProd_[0*8+4] = SymProd_[4*8+0] = 4;
    SymProd_[0*8+5] = SymProd_[5*8+0] = 5;
    SymProd_[0*8+6] = SymProd_[6*8+0] = 6;
    SymProd_[0*8+7] = SymProd_[7*8+0] = 7;
    SymProd_[1*8+2] = SymProd_[2*8+1] = 3;
    SymProd_[1*8+3] = SymProd_[3*8+1] = 2;
    SymProd_[1*8+4] = SymProd_[4*8+1] = 5;
    SymProd_[1*8+5] = SymProd_[5*8+1] = 4;
    SymProd_[1*8+6] = SymProd_[6*8+1] = 7;
    SymProd_[1*8+7] = SymProd_[7*8+1] = 6;
    SymProd_[2*8+3] = SymProd_[3*8+2] = 1;
    SymProd_[2*8+4] = SymProd_[4*8+2] = 6;
    SymProd_[2*8+5] = SymProd_[5*8+2] = 7;
    SymProd_[2*8+6] = SymProd_[6*8+2] = 4;
    SymProd_[2*8+7] = SymProd_[7*8+2] = 5;
    SymProd_[3*8+4] = SymProd_[4*8+3] = 7;
    SymProd_[3*8+5] = SymProd_[5*8+3] = 6;
    SymProd_[3*8+6] = SymProd_[6*8+3] = 5;
    SymProd_[3*8+7] = SymProd_[7*8+3] = 4;
    SymProd_[4*8+5] = SymProd_[5*8+4] = 1;
    SymProd_[4*8+6] = SymProd_[6*8+4] = 2;
    SymProd_[4*8+7] = SymProd_[7*8+4] = 3;
    SymProd_[5*8+6] = SymProd_[6*8+5] = 3;
    SymProd_[5*8+7] = SymProd_[7*8+5] = 2;
    SymProd_[6*8+7] = SymProd_[7*8+6] = 1;

    int * tmp_dim;
    tmp_dim = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)tmp_dim,'\0',nirrep_*sizeof(int));    

/* this is Q-Chem related 
    for ( int i=0; i < nmo_; i++) {

         orbital_symmetries_[i] = symmetry_space_energy_order_[i] - 1;

    }
*/

    // energy -> class & energy -> irrep & class -> irrep map

    int e2c_ind = 0;

    for (int h = 0; h < nirrep_; h++){

        for (int i=0;i<nrstc_;i++){

            if ( orbital_symmetries_[i] != h ) continue;

            OindMap_e2i_[i]       = tmp_dim[h];
            OindMap_e2c_[i]       = e2c_ind;
            OindMap_c2e_[e2c_ind] = i; 
            OindMap_c2i_[e2c_ind] = tmp_dim[h];

            tmp_dim[h]++;
            e2c_ind++;
               
        }

    }

    for (int h = 0; h < nirrep_; h++){

        for (int i=nrstc_; i < amo_ + nrstc_; i++){

            if ( orbital_symmetries_[i] != h ) continue;

            OindMap_e2i_[i] = tmp_dim[h];
            OindMap_e2c_[i] = e2c_ind;
            OindMap_c2e_[e2c_ind] = i;
            OindMap_c2i_[e2c_ind] = tmp_dim[h];

            tmp_dim[h]++;
            e2c_ind++;
               
        }

    }


    for (int h = 0;h < nirrep_; h++){

        for (int i=amo_+nrstc_; i < nrstv_ + amo_ + nrstc_; i++){

            if ( orbital_symmetries_[i] != h ) continue;

            OindMap_e2i_[i] = tmp_dim[h];
            OindMap_e2c_[i] = e2c_ind;
            OindMap_c2e_[e2c_ind] = i;
            OindMap_c2i_[e2c_ind] = tmp_dim[h];
              
            tmp_dim[h]++;
            e2c_ind++;

        }

    }

    // Geminal mappings

    memset((void*)tmp_dim,'\0',nirrep_*sizeof(int));

    for (int i = nrstc_; i < amo_ + nrstc_; i++){
 
        int h_i = orbital_symmetries_[i];

        for ( int j = nrstc_; j <= i; j++){

            int h_j = orbital_symmetries_[j];

            int h_ij = GemSym(h_i,h_j);

            Active_GindMap_e_[Full_ij_index(i,j,nmo_)] = tmp_dim[h_ij];
            Active_GindMap_e_[Full_ij_index(j,i,nmo_)] = tmp_dim[h_ij];

            tmp_dim[h_ij]++;

        }

    }

    // Class index matrix

    int nmo_have = 0;

    for ( int h = 0; h < nirrep_; h++ ){

        first_index_[h][0] = nmo_have;
        last_index_[h][0]  = first_index_[h][0]+rstcpi_[h];
        nmo_have += rstcpi_[h];

    }
 
    for ( int h = 0; h < nirrep_; h++ ){

        first_index_[h][1] = nmo_have;
        last_index_[h][1]  = first_index_[h][1]+amopi_[h];
        nmo_have += amopi_[h];

    }

    for ( int h = 0; h < nirrep_; h++ ){

        first_index_[h][2] = nmo_have;
        last_index_[h][2]  = first_index_[h][2]+rstvpi_[h];
        nmo_have += rstvpi_[h];

    }

    int p_df = 0;

    for ( int h_p = 0; h_p < nirrep_; h_p++){

        for ( int p_class = 0; p_class < 3; p_class++){

            for ( int p_c = first_index_[h_p][p_class]; p_c < last_index_[h_p][p_class]; p_c++){

                OindMap_c2df_[p_c] = p_df;

                p_df++;

            }
     
        }
   
    }

    for ( int h = 1; h < nirrep_; h++){

        full_oe_offset_[h] = full_oe_offset_[h-1] + nmopi_[h-1] * nmopi_[h-1];
        lt_oe_offset_[h]   = lt_oe_offset_[h-1] + nmopi_[h-1] * ( nmopi_[h-1] + 1 ) / 2;    

    }

    for ( int h_ij = 0; h_ij < nirrep_; h_ij++){

        int ngem_tmp = 0;

        for ( int h_i=0; h_i<nirrep_; h_i++){

            int h_j = GemSym(h_ij,h_i); 
            if ( h_j > h_i ) continue;
            
            for ( int i_c = first_index_[h_i][1]; i_c < last_index_[h_i][1]; i_c++){

                int j_c_max = last_index_[h_j][1];
                if ( h_i == h_j ) j_c_max = i_c + 1;

                for ( int j_c = first_index_[h_j][1]; j_c < j_c_max; j_c++){

                    Active_GindMap_c_[i_c * nmo_ + j_c] = ngem_tmp;
                    Active_GindMap_c_[j_c * nmo_ + i_c] = ngem_tmp;                     

                    ngem_tmp++;

                }
 
            }

        }

    }

    memset((void*)tmp_dim,'\0',nirrep_*sizeof(int));

    for ( int p_c = 0; p_c < nmo_; p_c++){

        int h_p = orbital_symmetries_[OindMap_c2e_[p_c]];

        for ( int q_c = 0; q_c <= p_c ; q_c++){

            int h_q = orbital_symmetries_[OindMap_c2e_[q_c]];
            int h_pq = GemSym(h_p,h_q);

            Full_GindMap_c_[p_c * nmo_ + q_c] = tmp_dim[h_pq];
            Full_GindMap_c_[q_c * nmo_ + p_c] = tmp_dim[h_pq];

            tmp_dim[h_pq]++;

        }

    }

    // number ofgeminals

    for ( int h = 0; h< nirrep_; h++ ){

        for ( int h_i = 0; h_i < nirrep_; h_i++){

            int h_j = GemSym(h,h_i);

            if ( h_j < h_i ){

                doc_gempi_[h] += rstcpi_[h_i]*rstcpi_[h_j];
                act_gempi_[h] += amopi_[h_i]*amopi_[h_j];
                ext_gempi_[h] += rstvpi_[h_i]*rstvpi_[h_j];

            }

            else if ( h_j == h_i ){

                doc_gempi_[h] += rstcpi_[h_i]*(rstcpi_[h_i]+1)/2;
                act_gempi_[h] += amopi_[h_i]*(amopi_[h_i]+1)/2;
                ext_gempi_[h] += rstvpi_[h_i]*(rstvpi_[h_i]+1)/2;

            }

        }

        
    }

    for ( int h_tu =1; h_tu < nirrep_; h_tu++) d2_irrep_offset_[h_tu] = d2_irrep_offset_[h_tu-1] +\
                                                act_gempi_[h_tu-1] * ( act_gempi_[h_tu-1] + 1 ) / 2;

    nnz_Q_ = amopi_[0]*nmopi_[0];

    for ( int h = 1; h < nirrep_; h++){

        Q_offset_[h] = Q_offset_[h-1] + amopi_[h-1]*nmopi_[h-1];
        
        nnz_Q_ += amopi_[h]*nmopi_[h];
    }

    for ( int p_e = 0; p_e < nmo_; p_e++){

        int h_p = orbital_symmetries_[p_e];

        for ( int q_e = 0; q_e <= p_e; q_e++){

             int h_q = orbital_symmetries_[q_e];
             int h_pq = GemSym(h_p,h_q);

             full_gempi_[h_pq]++;

        }

    }

    free(tmp_dim); 

}

int OrbitalOptimizer::GemSym( int h_i, int h_j){
    return SymProd_[h_i*8 + h_j]; 
}

int OrbitalOptimizer::Full_ij_index( int row, int col, int nrow){
    return col * nrow + row;
}

int OrbitalOptimizer::Lt_ij_index( int row, int col){

    int index;

    if ( row > col) {

       index = row * ( row + 1 )/2 + col;        

    }
    else {   
 
       index = col * ( col + 1 )/2 + row;
 
    }

    return index;

}

int OrbitalOptimizer::Lt_ij_oe_index( int p_c, int q_c, int h){

    int p_i = OindMap_c2i_[p_c];
    int q_i = OindMap_c2i_[q_c];
   
    int index = Lt_ij_index(p_i, q_i) + lt_oe_offset_[h];

    return index;

}

int OrbitalOptimizer::Full_ij_oe_index( int row_c, int col_c, int h){

    int row_i = OindMap_c2i_[row_c];
    int col_i = OindMap_c2i_[col_c];

    int index = Full_ij_index(row_i, col_i,nmo_) + full_oe_offset_[h];

    return index;

}

int OrbitalOptimizer::df_ij_index( int p_c, int q_c ){

    int p_df = OindMap_c2df_[p_c];
    int q_df = OindMap_c2df_[q_c];

    int index = Lt_ij_index(p_df,q_df);

    return index;

}

int OrbitalOptimizer::Q_index( int t_c, int p_c, int h_p ){

    int t_red = t_c - first_index_[h_p][1]; 
    int p_i   = OindMap_c2i_[p_c] ;

    int index = Q_offset_[h_p] + p_i * amopi_[h_p] + t_red;

    return index;

}



int OrbitalOptimizer::Max_value( int * array, int dim){

    int maxval = 0;

    for ( int i = 0; i < dim; i++)

        if ( maxval < array[i] )
 
            maxval = array[i];

    return maxval;
}



void OrbitalOptimizer::Alloc_OindMap(){

    int nmo = nrstc_ + amo_ + nrstv_;

    SymProd_ = (int*)malloc(64*sizeof(int));
    memset((void*)SymProd_,'\0',64*sizeof(int));

    doc_gempi_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)doc_gempi_,'\0',nirrep_*sizeof(int));

    act_gempi_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)act_gempi_,'\0',nirrep_*sizeof(int));       

    ext_gempi_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)ext_gempi_,'\0',nirrep_*sizeof(int));

    full_gempi_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)full_gempi_,'\0',nirrep_*sizeof(int));

    lt_oe_offset_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)lt_oe_offset_,'\0',nirrep_*sizeof(int));

    full_oe_offset_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)full_oe_offset_,'\0',nirrep_*sizeof(int));    

    d2_irrep_offset_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)d2_irrep_offset_,'\0',nirrep_*sizeof(int));

    Q_offset_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)Q_offset_,'\0',nirrep_*sizeof(int));

    aa_ras1pi_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)aa_ras1pi_,'\0',nirrep_*sizeof(int));

    OindMap_e2c_ = (int *)malloc(nmo*sizeof(int)); 
    memset((void*)OindMap_e2c_,'\0',nmo*sizeof(int));

    OindMap_c2e_ = (int *)malloc(nmo*sizeof(int));
    memset((void*)OindMap_c2e_,'\0',nmo*sizeof(int));

    OindMap_e2i_ = (int *)malloc(nmo*sizeof(int));
    memset((void*)OindMap_e2i_,'\0',nmo*sizeof(int));

    OindMap_c2i_ = (int *)malloc(nmo*sizeof(int));
    memset((void*)OindMap_c2i_,'\0',nmo*sizeof(int));

    OindMap_c2df_ = (int *)malloc(nmo*sizeof(int));
    memset((void*)OindMap_c2df_,'\0',nmo*sizeof(int));

    Active_GindMap_e_ = (int *)malloc(nmo*nmo*sizeof(int));
    memset((void*)Active_GindMap_e_,'\0',nmo*nmo*sizeof(int));

    Active_GindMap_c_ = (int *)malloc(nmo*nmo*sizeof(int));
    memset((void*)Active_GindMap_c_,'\0',nmo*nmo*sizeof(int));

    Full_GindMap_c_ = (int *)malloc(nmo*nmo*sizeof(int));
    memset((void*)Full_GindMap_c_,'\0',nmo*nmo*sizeof(int));

    first_index_ = (int**)malloc(nirrep_*sizeof(int*));
    for (int i = 0; i < nirrep_; i++){
        first_index_[i] = (int*)malloc(3*sizeof(int));
    }

    last_index_ = (int**)malloc(nirrep_*sizeof(int*));
    for (int i = 0; i < nirrep_; i++){
        last_index_[i] = (int*)malloc(3*sizeof(int));
    }
};

void OrbitalOptimizer::Dealloc_OindMap(){

    free(SymProd_);
    free(doc_gempi_);
    free(act_gempi_);
    free(ext_gempi_);
    free(d2_irrep_offset_);
    free(full_gempi_);
    free(full_oe_offset_);
    free(Q_offset_);
    free(lt_oe_offset_);
    free(OindMap_e2c_);           
    free(OindMap_c2e_);
    free(OindMap_e2i_);
    free(OindMap_c2i_);
    free(OindMap_c2df_);
    free(Active_GindMap_c_);
    free(Active_GindMap_e_);
    free(Full_GindMap_c_);
    free(aa_ras1pi_); 

    for (int i = 0; i < nirrep_; i++){
        free(first_index_[i]);
    }
    free(first_index_);

    for (int i = 0; i < nirrep_; i++){
        free(last_index_[i]);
    }
    free(last_index_);


};

} // end of namespace
