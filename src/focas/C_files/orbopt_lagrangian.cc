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

void OrbitalOptimizer::Get_GeneralizedFock(double * d2, double * d1, double * tei, double * oei, double * Fock){

    // function to calculate generalized Fock matrix

/*
    // these functions are now called in init

    // Allocate/initialize indexing arrays
    Alloc_OindMap();
  
    Determine_OindMap();
*/
 
    // sort oei, d1, and d2 (needs maping arrays)
    int direction = 1;
    Sort( d2, d1, oei, direction );
    
    // Allocate gradient arrays
    Alloc_Gradient();    

    // Calculate the generalized Fock matrix
    Compute_GeneralizedFock(d2, d1, tei, oei, Fock);

    // Unsort d2, d1, and oei arrays
    direction = -1;
    Sort( d2, d1, oei, direction );

    // Deallocate arrays
    Dealloc_Gradient();

/*

    // this function is now called in destructor
    Dealloc_OindMap();

*/

}

void OrbitalOptimizer::Compute_GeneralizedFock(double * d2, double * d1, double * tei, double * oei, double * Fock){

    Initialize_scratch_G();

    // Compute inactive Fock matrix
    Compute_Fi( tei, oei );
 
    // Compute active Fock matrix
    Compute_Fa( tei, d1 );

    // Compute contractions of d2 and tei
    Compute_Q( tei, d2 );
 
    // compute contraction of d1 and Fi
    Compute_Z( d1 );

    for ( int p_class = 0; p_class < 2; p_class++ ) {

        for ( int q_class = 0; q_class < 3; q_class++ ){

            GeneralizedFock_block(Fock, p_class,q_class);

        }   

    }

}

void OrbitalOptimizer::GeneralizedFock_block(double * Fock, int p_class, int q_class){

    int offset = 0;
  
    for ( int p_class = 0; p_class < 2; p_class++ ) {

        for ( int q_class = 0; q_class < 3; q_class++ ){
    
            for ( int h = 0; h < nirrep_; h++ ){

                int nmo_q = last_index_[h][q_class] - first_index_[h][q_class];
                int nmo_p = last_index_[h][p_class] - first_index_[h][p_class];

                if ( nmo_p * nmo_q == 0 ) continue;

                int p = 0;

                for ( int p_c = first_index_[h][p_class]; p_c < last_index_[h][p_class]; p_c++ ){

                    int q = 0;

                    for ( int q_c = first_index_[h][q_class]; q_c < last_index_[h][q_class]; q_c++ ){

                        double Fval=0.0e0;

                        if ( p_class == 0 ) {

                            // d-d, a-d, and e-d 
  
                            Fval = 2.0e0 * ( Fi_[Lt_ij_oe_index(q_c,p_c,h) ] + \
                                             Fa_[Lt_ij_oe_index(q_c,p_c,h) ] );

                        }

                        else if ( p_class == 1 ){

                            // a-a and e-a
                             
                            Fval = Q_[ Q_index(p_c,q_c,h) ] + Z_[ Q_index(p_c,q_c,h) ];
                     
                        }

                        Fock[ offset + p * nmo_q + q ] = Fval;
                        
                        q++;

                    }

                    p++;

                }

                offset += nmo_q * nmo_p;

            }

        }

    }

}

}
