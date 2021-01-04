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

void OrbitalOptimizer::ExponentiateKappa(double * kappa){

    for ( int h = 0; h < nirrep_; h++){

        if ( nmopi_[h] == 0 ) continue;

        ExponentiateKappa_UnpackSymBlock(kappa, h);

        ExponentiateKappa_ComputeUSymBlock(h);

    }

}

void OrbitalOptimizer::ExponentiateKappa_ComputeUSymBlock(int h){

    int nmo    = nmopi_[h];

    int U_off  = U_offset_[h];

    int lwork  = (1 + 6 * nmo + 2 * nmo * nmo);

    int liwork = 3 + 5 * nmo;

    // zero this block of U

    memset((void*)&U_[U_off],'\0',nmo*nmo*sizeof(double));

    // allocate tmeporary matrices

    double * K2 = (double *)malloc(nmo*nmo*sizeof(double));
    memset((void*)K2,'\0',nmo*nmo*sizeof(double));

    double * d = (double *)malloc(nmo*nmo*sizeof(double));
    memset((void*)d,'\0',nmo*nmo*sizeof(double));

    double * tmp_1 = (double *)malloc(nmo*nmo*sizeof(double));
    memset((void*)tmp_1,'\0',nmo*nmo*sizeof(double));

    double * tmp_2 = (double *)malloc(nmo*nmo*sizeof(double));
    memset((void*)tmp_2,'\0',nmo*nmo*sizeof(double));

    double * work = (double *)malloc(lwork*sizeof(double));
    memset((void*)work,'\0',lwork*sizeof(double));

    int * iwork   = (int*)malloc(liwork * sizeof(int));
    memset((void*)iwork,'\0',liwork*sizeof(int));

    // calcuate Kappa**2 and store in K2

    F_DGEMM('n','n',nmo,nmo,nmo,1.0e0,&Tei_scr1_[0],nmo, \
                                      &Tei_scr1_[0],nmo, \
                                0.0e0,&K2[0],nmo);

    // diagonalize Kappa**2

    C_DSYEVD('v','u',nmo,K2,nmo,d,work,lwork,iwork,liwork);

    // ensure that eigenvalues are positive

    for ( int p = 0; p < nmo; p++){

        double val = d[p];

        if ( val < 0.0e0 ) val = -val;

        d[p] = sqrt(val);

    }

    // Compute and store X * sin(d) * d^-1 in tmp_1

    int offset = 0;
    for ( int p = 0; p < nmo; p++ ){

        double val = 1.0e0;
        if ( d[p] != 0.0e0 ) val = sin(d[p]) / d[p];

        C_DCOPY(nmo,&K2[offset],1,&tmp_1[offset],1);
        C_DSCAL(nmo,val,&tmp_1[offset],1);

        offset += nmo;

    }

    // Compute and store X * d^-1 * sin(d) X^T in tmp_2

    F_DGEMM('n','t',nmo,nmo,nmo,1.0e0,&K2[0],nmo,     \
                                      &tmp_1[0],nmo, \
                                0.0e0,&tmp_2[0],nmo);

    // Compute and store K * X * d^-1 * sin(d) X^T in U_

    F_DGEMM('n','n',nmo,nmo,nmo,1.0e0,&Tei_scr1_[0],nmo, \
                                      &tmp_2[0],nmo,     \
                                0.0e0,&U_[U_off],nmo);

    // compute X*cos(d)

    offset = 0;
    for ( int p = 0; p < nmo; p++){

        double val = cos(d[p]);

        C_DCOPY(nmo,&K2[offset],1,&tmp_1[offset],1);
        C_DSCAL(nmo,val,&tmp_1[offset],1);

        offset += nmo;

    }

    // compute and add X * cos(d) * X^T to U_

    F_DGEMM('n','t',nmo,nmo,nmo,1.0e0,&K2[0],nmo,      \
                                      &tmp_1[0],nmo,  \
                                1.0e0,&U_[U_off],nmo);

    free(tmp_2);
    free(tmp_1);
    free(iwork);
    free(work);
    free(d);
    free(K2);

}

void OrbitalOptimizer::ExponentiateKappa_UnpackSymBlock(double * kappa, int h){

    memset((void*)Tei_scr1_,'\0',Tei_scr_dim_*sizeof(double));

    int nmo_p = nmopi_[h];

    for ( int p_class = 0; p_class < 3; p_class ++ ){

        int q_class_min = p_class + 1;
        if ( active_active_rotations_ && p_class == 1 ) q_class_min = p_class;

        for ( int q_class = q_class_min; q_class < 3; q_class++){

            for ( int p_c = first_index_[h][p_class]; p_c < last_index_[h][p_class]; p_c++){

                int p_i = OindMap_c2i_[p_c];

                int q_min = first_index_[h][q_class];
                if ( q_class == p_class ) q_min = p_c + 1;

                for ( int q_c = q_min; q_c < last_index_[h][q_class]; q_c++){

                    int q_i = OindMap_c2i_[q_c];

                    int grad_ind = Grad_IndMap_c_[Lt_ij_index(p_c,q_c)];

                    Tei_scr1_[Full_ij_index(q_i,p_i,nmo_p)] =  kappa[grad_ind];
                    Tei_scr1_[Full_ij_index(p_i,q_i,nmo_p)] = -kappa[grad_ind];

                }

            }

        }

    }
}

void OrbitalOptimizer::Alloc_U(){

    U_offset_ = (int *)malloc(nirrep_*sizeof(int));
    memset((void*)U_offset_,'\0',nirrep_*sizeof(int));

    int U_dim_    = 0;
    Tei_scr_dim_ = 0;

    for ( int h = 0; h < nirrep_; h++){

        U_offset_[h] = U_dim_;

        U_dim_ += nmopi_[h]*nmopi_[h];

        if ( Tei_scr_dim_ < nmopi_[h]*nmopi_[h] ) Tei_scr_dim_ = nmopi_[h]*nmopi_[h];

    }

    U_ = (double *)malloc(U_dim_*sizeof(double));
    memset((void*)U_,'\0',U_dim_*sizeof(double));

    Tei_scr1_ = (double *)malloc(max_thread_*Tei_scr_dim_*sizeof(double));
    memset((void*)Tei_scr1_,'\0',max_thread_*Tei_scr_dim_*sizeof(double));

    Tei_scr2_ = (double *)malloc(max_thread_*Tei_scr_dim_*sizeof(double));
    memset((void*)Tei_scr2_,'\0',max_thread_*Tei_scr_dim_*sizeof(double));

}

void OrbitalOptimizer::Dealloc_U(){

    free(U_offset_);
    free(U_);
    free(Tei_scr1_);
    free(Tei_scr2_);

}



}
