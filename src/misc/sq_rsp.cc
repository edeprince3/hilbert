/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
** \file
** \brief Diagnoalize a symmetrix square matrix
** \ingroup CIOMR
*/

#include "psi4/psifiles.h"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "sq_rsp.h"

#include <cstdlib>
#include <cstring>

using namespace psi;

namespace hilbert {

/* translation into c of a translation into FORTRAN77 of the EISPACK */
/* matrix diagonalization routines */

/**
*  WARNING: Psi 3 Fortran routine sq_rsp  deprecated
*  by Robert Parrish, robparrish@gmail.com
*
*  sq_rsp now calls the LAPACK method DSYEV for
*  numerical stability, speed, and threading
*
*  the signature of this method remains the same
*
*  June 22, 2010
**/

/*!
** sq_rsp(): diagomalize a symmetric square matrix ('array').
**
** \param nm     = rows of matrix
** \param n      = columns of matrix
** \param array  = matrix to diagonalize
** \param e_vals = array to hold eigenvalues
** \param matz   = 0 (no eigenvectors, eigenvals in ascending order)
**               = 1 (eigenvectors and eigenvalues in ascending order)
**               = 2 (no eigenvectors, eigenvalues in descending order)
**               = 3 (eigenvectors and eigenvalues in descending order)
** \param e_vecs = matrix of eigenvectors (one column for each eigvector)
** \param toler  = tolerance for eigenvalues?  Often 1.0E-14.
**
** \ingroup CIOMR
*/
void sq_rsp(int /*nm*/, int n, double** array, double* e_vals, int matz, double** e_vecs, double /*toler*/) {
    if ((matz > 3) || (matz < 0)) matz = 0;
    // Do you want eigenvectors?
    bool eigenvectors = (matz == 1 || matz == 3) ? true : false;
    // Ascending or Descending?
    bool ascending = (matz == 0 || matz == 1) ? true : false;

    // Get Eigenvectors (Use 'V')
    if (eigenvectors) {
        // First Temp array (the lda isn't so fly in libciomr)
        double** Temp_sqrsp = block_matrix(n, n);

        // Copy array to Temp_sqrsp (for loops required)
        for (int r = 0; r < n; r++)
            for (int c = 0; c < n; c++) Temp_sqrsp[r][c] = array[r][c];

        // outfile->Printf("  Initial Matrix:\n");
        // print_mat(e_vecs,n,n,outfile);
        // printf("%d",n);

        // Get scratch vector and call DSYEV
        // The eigenvectors are placed in e_vecs in ascending order
        int lwork_sqrsp  = (1 + 6 * n + 2 * n * n);
        double* work_sqrsp = (double*)malloc(lwork_sqrsp*sizeof(double));

        int liwork = 3 + 5 * n;
        int * iwork   = (int*)malloc(liwork * sizeof(int));
        memset((void*)iwork,'\0',liwork*sizeof(int));

        C_DSYEVD('V', 'U', n, &Temp_sqrsp[0][0], n, &e_vals[0], &work_sqrsp[0], lwork_sqrsp, iwork, liwork);
        free(work_sqrsp);
        free(iwork);


        // outfile->Printf("  Eigenvectors:\n");
        // print_mat(e_vecs,n,n,outfile);

        // printf("\n  Eigenvalues:\n");
        // for (int r = 0; r<n; r++)
        //    printf("  r = %d, %14.10f\n",r+1,e_vals[r]);

        // LAPACK stores eigenvectors in rows, we need them in columns
        // This overhead is why you should always call DSYEV!
        double** T_sqrsp_2 = block_matrix(n, n);
        C_DCOPY(static_cast<size_t>(n) * n, &Temp_sqrsp[0][0], 1, &T_sqrsp_2[0][0], 1);
        for (int r = 0; r < n; r++) C_DCOPY(n, &T_sqrsp_2[r][0], 1, &Temp_sqrsp[0][r], n);
        free_block(T_sqrsp_2);

        // If descending is required, the canonical order must be reversed
        // Sort is stable
        if (!ascending) {
            double* Temp_sqrsp_col = (double*)malloc(n*sizeof(double));
            double w_Temp_sqrsp;

            for (int c = 0; c < n / 2; c++) {
                // Swap eigenvectors
                C_DCOPY(n, &Temp_sqrsp[0][c], n, &Temp_sqrsp_col[0], 1);
                C_DCOPY(n, &Temp_sqrsp[0][n - c - 1], n, &Temp_sqrsp[0][c], n);
                C_DCOPY(n, &Temp_sqrsp_col[0], 1, &Temp_sqrsp[0][n - c - 1], n);

                // Swap eigenvalues
                w_Temp_sqrsp = e_vals[c];
                e_vals[c] = e_vals[n - c - 1];
                e_vals[n - c - 1] = w_Temp_sqrsp;
            }

            free(Temp_sqrsp_col);
        }
        // Copy from Temp_sqrsp to e_vecs (for loops required)
        for (int r = 0; r < n; r++)
            for (int c = 0; c < n; c++) e_vecs[r][c] = Temp_sqrsp[r][c];

        free_block(Temp_sqrsp);
        // outfile->Printf("  Eigenvectors (After sort):\n");
        // print_mat(e_vecs,n,n,outfile);

        // printf("\n  Eigenvalues (After sort):\n");
        // for (int r = 0; r<n; r++)
        //    printf("  r = %d, %14.10f\n",r+1,e_vals[r]);

        // No Eigenvectors (Use 'N')
    } else {
        // e_vecs might not be initialized
        // Use a Temp array
        double** Temp_sqrsp = block_matrix(n, n);
        // Copy array to Temp
        // must use for loops because of bloody irrep blocking in CSCF
        for (int r = 0; r < n; r++)
            for (int c = 0; c < n; c++) Temp_sqrsp[r][c] = array[r][c];

        // Form scratch array and call DSYEV
        //'N' in parameter 1 to DSYEV terminates
        // the algorithm after eigenvectors are found
        // Canonical order is ascending in LAPACK
        int lwork_sqrsp  = (1 + 6 * n + 2 * n * n);
        double* work_sqrsp = (double*)malloc(lwork_sqrsp*sizeof(double));

        int liwork = 3 + 5 * n;
        int * iwork   = (int*)malloc(liwork * sizeof(int));
        memset((void*)iwork,'\0',liwork*sizeof(int));

        C_DSYEVD('N', 'U', n, &Temp_sqrsp[0][0], n, &e_vals[0], &work_sqrsp[0], lwork_sqrsp, iwork, liwork);
        free(iwork);
        free(work_sqrsp);
        free_block(Temp_sqrsp);

        // If descending is required, the canonical order must be reversed
        // Sort is stable
        if (!ascending) {
            double w_Temp_sqrsp;

            // Eigenvalues only
            for (int c = 0; c < n / 2; c++) {
                w_Temp_sqrsp = e_vals[c];
                e_vals[c] = e_vals[n - c - 1];
                e_vals[n - c - 1] = w_Temp_sqrsp;
            }
        }
    }
}
}
