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

#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/mintshelper.h>

#include <psi4/libtrans/integraltransform.h>
#include <psi4/libtrans/mospace.h>

#include <psi4/libmints/writer.h>
#include <psi4/libmints/writer_file_prefix.h>
#include <psi4/libmints/molecule.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/vector.h>
#include <psi4/libmints/matrix.h>

#include "v2rdm_solver.h"

#include <misc/omp.h>

using namespace psi;

namespace hilbert{


// 
// dump 1-/2-electron integrals and 1-/2-RDM to disk
// 
// notes:
//
// - all quantities should be in the NO basis
// - all 1-/2-electron integrals are required
// - only active 1-/2-RDM elements are required
// - as a first pass, this will only work with DF integrals
//
void v2RDMSolver::FCIDUMP() {

    if ( nfrzv_ > 0 ) {
        throw PsiException("FCIDUMP does not work with frozen virtual orbitals",__FILE__,__LINE__);
    }

    std::string filename = get_writer_file_prefix(reference_wavefunction_->molecule()->name());
    std::string rdm_filename = filename + ".rdm";
    std::string int_filename = filename + ".int";

    FILE * int_fp = fopen(int_filename.c_str(),"wb");
    FILE * rdm_fp = fopen(rdm_filename.c_str(),"wb");

    int zero = 0;

    // count number of inactive orbitals:
    int ninact = 0;
    for (int h = 0; h < nirrep_; h++) {
        ninact += frzcpi_[h];
        ninact += rstcpi_[h];
    }

    // orbitals are ordered by irrep and by space within each irrep. 
    // need a map to order by space and by irrep within each space
    // (plus 1)
    int * map = (int*)malloc(nmo_*sizeof(int));
    int count = 0;
    // core
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h] + rstcpi_[h]; i++) {
            map[i + pitzer_offset_full[h]] = count + 1;
            count++;
        }
    }
    // active
    for (int h = 0; h < nirrep_; h++) {
        for (int t = 0; t < amopi_[h]; t++) {
            map[t + pitzer_offset_full[h] + frzcpi_[h] + rstcpi_[h]] = count + 1;
            count++;
        }
    }
    // virtual
    for (int h = 0; h < nirrep_; h++) {
        for (int a = 0; a < nmopi_[h] - ( frzcpi_[h] + rstcpi_[h] + amopi_[h] ); a++) {
            map[a + pitzer_offset_full[h] + frzcpi_[h] + rstcpi_[h] + amopi_[h]] = count + 1;
            count++;
        }
    }

    // two-electron integrals
    for (int p = 0; p < nmo_; p++) {
        int hp = symmetry_really_full[p];
        for (int q = p; q < nmo_; q++) {
            int hq = symmetry_really_full[q];
            int hpq = hp ^ hq;
            long int pq = INDEX(p,q);
            for (int r = 0; r < nmo_; r++) {
                int hr = symmetry_really_full[r];
                for (int s = r; s < nmo_; s++) {
                    int hs = symmetry_really_full[s];
                    int hrs = hr ^ hs;
                    if ( hpq != hrs ) continue;
                    long int rs = INDEX(r,s);
                    if ( pq > rs ) continue;
                    double dum = TEI(p,q,r,s,hpq);
                    //if ( fabs(dum) < 1e-12 ) continue;
                    int pp = map[p];
                    int qq = map[q];
                    int rr = map[r];
                    int ss = map[s];
                    fwrite (&dum , sizeof(double), 1, int_fp);
                    fwrite (&pp , sizeof(int), 1, int_fp);
                    fwrite (&qq , sizeof(int), 1, int_fp);
                    fwrite (&rr , sizeof(int), 1, int_fp);
                    fwrite (&ss , sizeof(int), 1, int_fp);
                }
            }
        }
    }
    // one-electron integrals

    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));
    std::shared_ptr<Matrix> T (new Matrix(mints->so_kinetic()));
    std::shared_ptr<Matrix> V (new Matrix(mints->so_potential()));

    T->transform(Ca_);
    V->transform(Ca_);

    for (int h = 0; h < nirrep_; h++) {
        double ** Tp = T->pointer(h);
        double ** Vp = V->pointer(h);
        for (int p = 0; p < nmopi_[h]; p++) {
            for (int q = 0; q < nmopi_[h]; q++) {
                double dum = Tp[p][q] + Vp[p][q];
                //if ( fabs(dum) < 1e-12 ) continue;
                int pp = map[p + pitzer_offset_full[h]];
                int qq = map[q + pitzer_offset_full[h]];
                fwrite (&dum , sizeof(double), 1, int_fp);
                fwrite (&pp , sizeof(int), 1, int_fp);
                fwrite (&qq , sizeof(int), 1, int_fp);
                fwrite (&zero , sizeof(int), 1, int_fp);
                fwrite (&zero , sizeof(int), 1, int_fp);
            }
        }
    }
    fwrite(&enuc_, sizeof(double), 1, int_fp);
    fwrite(&zero, sizeof(int), 1, int_fp);
    fwrite(&zero, sizeof(int), 1, int_fp);
    fwrite(&zero, sizeof(int), 1, int_fp);
    fwrite(&zero, sizeof(int), 1, int_fp);
    fclose(int_fp);

    // two-electron rdm
    double * x_p = x->pointer();
    double e2 = 0.0;
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];

            int ji = ibas_ab_sym[h][j][i];

            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];

                int lk = ibas_ab_sym[h][l][k];

                double dum = x_p[d2aboff[h] + ij * gems_ab[h] + kl]
                           + x_p[d2aboff[h] + ji * gems_ab[h] + lk];

                if ( i != j && k != l ) {

                    int ija = ibas_aa_sym[h][i][j];
                    int kla = ibas_aa_sym[h][k][l];
                    int sg = 1;
                    if ( i > j ) sg = -sg;
                    if ( k > l ) sg = -sg;

                    dum += sg * x_p[d2aaoff[h] + ija * gems_aa[h] + kla];
                    dum += sg * x_p[d2bboff[h] + ija * gems_aa[h] + kla];

                }

                int hi  = symmetry[i];
                int hk  = symmetry[k];
                int hik = hi ^ hk;
                e2 += 0.5 * dum * TEI(full_basis[i],full_basis[k],full_basis[j],full_basis[l],hik);

                //if ( fabs(dum) < 1e-12 ) continue;
                int ii = map[full_basis[i]] - ninact;
                int jj = map[full_basis[j]] - ninact;
                int kk = map[full_basis[k]] - ninact;
                int ll = map[full_basis[l]] - ninact;
                fwrite (&dum , sizeof(double), 1, rdm_fp);
                fwrite (&ii , sizeof(int), 1, rdm_fp);
                fwrite (&jj , sizeof(int), 1, rdm_fp);
                fwrite (&kk , sizeof(int), 1, rdm_fp);
                fwrite (&ll , sizeof(int), 1, rdm_fp);
            }
        }
    }

    double e1 = 0.0;

    // one-electron rdm
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];

            int hi = symmetry[i];
            int hj = symmetry[j];

            if ( hi != hj ) continue;

            int ii = i - pitzer_offset[hi];
            int jj = j - pitzer_offset[hj];

            double dum = x_p[d1aoff[hi] + ii * amopi_[hi] + jj] + x_p[d1boff[hi] + ii * amopi_[hi] + jj];

            //double ** Tp = T->pointer(hi);
            //double ** Vp = V->pointer(hi);

            //ii += rstcpi_[hi];
            //jj += rstcpi_[hi];
            //e1 += dum * (Tp[ii][jj] + Vp[ii][jj]);

            ii = map[full_basis[i]] - ninact;
            jj = map[full_basis[j]] - ninact;

            fwrite (&dum , sizeof(double), 1, rdm_fp);
            fwrite (&ii , sizeof(int), 1, rdm_fp);
            fwrite (&jj , sizeof(int), 1, rdm_fp);
            fwrite (&zero , sizeof(int), 1, rdm_fp);
            fwrite (&zero , sizeof(int), 1, rdm_fp);

        }
    }

    //for (int h = 0; h < nirrep_; h++) {
    //    double ** Tp = T->pointer(h);
    //    double ** Vp = V->pointer(h);
    //    for (int i = 0; i < rstcpi_[h] + frzcpi_[h]; i++) {
    //        e1 += 2.0 * (Tp[i][i] + Vp[i][i]);
    //    }
    //}
    //printf("%20.12lf\n",e1);

    fclose(rdm_fp);


}


double rdm_entropy(int nirrep, int * dim, double * x, int * off){

    SharedMatrix mat (new Matrix(nirrep,dim,dim));
    SharedMatrix eigvec (new Matrix(nirrep,dim,dim));
    SharedVector eigval (new Vector(nirrep,dim));

    double entropy = 0.0;

    for (int h = 0; h < nirrep; h++) {
        for (int i = 0; i < dim[h]; i++) {
            for (int j = 0; j < dim[h]; j++) {
                mat->pointer(h)[i][j]  = x[off[h]+i*dim[h]+j];
            }
        }
    }
    mat->diagonalize(eigvec,eigval,descending);

    for (int h = 0; h < nirrep; h++) {
        for (int i = 0; i < dim[h]; i++) {
            double dum = eigval->pointer(h)[i];
            if ( fabs(dum) < 1e-12 && dum < 0.0 ) dum = -dum;
            entropy -= dum * log2(dum);
        }
    }

    return entropy;
}

void v2RDMSolver::print_rdms() {

    std::string filename = get_writer_file_prefix(reference_wavefunction_->molecule()->name());
    std::string rdm_filename = filename + ".rdm";

    FILE * rdm_fp = fopen(rdm_filename.c_str(),"w");

    double * x_p = x->pointer();

    double na = nalpha_ - nrstc_ - nfrzc_;
    double nb = nbeta_ - nrstc_ - nfrzc_;

    int nact = 0;
    for (int h = 0; h < nirrep_; h++) {
        nact += amopi_[h];
    }
    double ms = (multiplicity_ - 1.0)/2.0;

    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    v2RDM @Nalpha           %20.12lf\n",na);
    fprintf(rdm_fp,"    v2RDM @Nbeta            %20.12lf\n",nb);
    fprintf(rdm_fp,"    v2RDM @Nact             %20.12lf\n",(double)nact);
    fprintf(rdm_fp,"    v2RDM @S2               %20.12lf\n",ms*(ms+1));
    fprintf(rdm_fp,"\n");

    // generalized entropy

    // D1

    fprintf(rdm_fp,"    v2RDM @entropy(D1a)     %20.12lf\n",rdm_entropy(nirrep_, amopi_, x_p, d1aoff));
    fprintf(rdm_fp,"    v2RDM @entropy(D1b)     %20.12lf\n",rdm_entropy(nirrep_, amopi_, x_p, d1boff));

    // Q1

    fprintf(rdm_fp,"    v2RDM @entropy(Q1a)     %20.12lf\n",rdm_entropy(nirrep_, amopi_, x_p, q1aoff));
    fprintf(rdm_fp,"    v2RDM @entropy(Q1b)     %20.12lf\n",rdm_entropy(nirrep_, amopi_, x_p, q1boff));

    // D2

    fprintf(rdm_fp,"    v2RDM @entropy(D2aa)    %20.12lf\n",rdm_entropy(nirrep_, gems_aa, x_p, d2aaoff));
    fprintf(rdm_fp,"    v2RDM @entropy(D2bb)    %20.12lf\n",rdm_entropy(nirrep_, gems_aa, x_p, d2bboff));
    fprintf(rdm_fp,"    v2RDM @entropy(D2ab)    %20.12lf\n",rdm_entropy(nirrep_, gems_ab, x_p, d2aboff));

    if ( constrain_q2_ ) {

        // Q2

        fprintf(rdm_fp,"    v2RDM @entropy(Q2aa)    %20.12lf\n",rdm_entropy(nirrep_, gems_aa, x_p, q2aaoff));
        fprintf(rdm_fp,"    v2RDM @entropy(Q2bb)    %20.12lf\n",rdm_entropy(nirrep_, gems_aa, x_p, q2bboff));
        fprintf(rdm_fp,"    v2RDM @entropy(Q2ab)    %20.12lf\n",rdm_entropy(nirrep_, gems_ab, x_p, q2aboff));

    }

    if ( constrain_g2_ ) {

        // G2

        int * tmp = (int*)malloc(nirrep_*sizeof(int));
        for (int h = 0; h < nirrep_; h++) {
            tmp[h] = 2 * gems_ab[h];
        }

        fprintf(rdm_fp,"    v2RDM @entropy(G2ab)    %20.12lf\n",rdm_entropy(nirrep_, gems_ab, x_p, g2aboff));
        fprintf(rdm_fp,"    v2RDM @entropy(G2ba)    %20.12lf\n",rdm_entropy(nirrep_, gems_ab, x_p, g2baoff));
        fprintf(rdm_fp,"    v2RDM @entropy(G2aa/bb) %20.12lf\n",rdm_entropy(nirrep_,     tmp, x_p, g2aaoff));

        free(tmp);

    }

    fprintf(rdm_fp,"\n");

    // norm of cumulant 2RDM

    std::shared_ptr<Matrix> del2aa (new Matrix(nirrep_, gems_aa, gems_aa));
    std::shared_ptr<Matrix> del2bb (new Matrix(nirrep_, gems_aa, gems_aa));
    std::shared_ptr<Matrix> del2ab (new Matrix(nirrep_, gems_ab, gems_ab));

    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            int hi = symmetry[i];
            int hj = symmetry[j];
            int ii = i - pitzer_offset[hi];
            int jj = j - pitzer_offset[hj];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                int hk = symmetry[k];
                int hl = symmetry[l];
                int kk = k - pitzer_offset[hk];
                int ll = l - pitzer_offset[hl];
                double dik = 0.0;
                double djl = 0.0;
                double dil = 0.0;
                double djk = 0.0;
                if ( hi == hk ) {
                    dik = x->pointer()[d1aoff[hi] + ii * amopi_[hi] + kk];
                }
                if ( hi == hl ) {
                    dil = x->pointer()[d1aoff[hi] + ii * amopi_[hi] + ll];
                }
                if ( hj == hk ) {
                    djk = x->pointer()[d1aoff[hj] + jj * amopi_[hj] + kk];
                }
                if ( hj == hl ) {
                    djl = x->pointer()[d1aoff[hj] + jj * amopi_[hj] + ll];
                }
                del2aa->pointer(h)[ij][kl] = x->pointer()[d2aaoff[h] + ij * gems_aa[h] + kl] - ( dik * djl - dil * djk );

            }
        }
    }

    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            int hi = symmetry[i];
            int hj = symmetry[j];
            int ii = i - pitzer_offset[hi];
            int jj = j - pitzer_offset[hj];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                int hk = symmetry[k];
                int hl = symmetry[l];
                int kk = k - pitzer_offset[hk];
                int ll = l - pitzer_offset[hl];
                double dik = 0.0;
                double djl = 0.0;
                double dil = 0.0;
                double djk = 0.0;
                if ( hi == hk ) {
                    dik = x->pointer()[d1boff[hi] + ii * amopi_[hi] + kk];
                }
                if ( hi == hl ) {
                    dil = x->pointer()[d1boff[hi] + ii * amopi_[hi] + ll];
                }
                if ( hj == hk ) {
                    djk = x->pointer()[d1boff[hj] + jj * amopi_[hj] + kk];
                }
                if ( hj == hl ) {
                    djl = x->pointer()[d1boff[hj] + jj * amopi_[hj] + ll];
                }
                del2bb->pointer(h)[ij][kl] = x->pointer()[d2bboff[h] + ij * gems_aa[h] + kl] - ( dik * djl - dil * djk );

            }
        }
    }

    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];
            int hi = symmetry[i];
            int hj = symmetry[j];
            int ii = i - pitzer_offset[hi];
            int jj = j - pitzer_offset[hj];
            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];
                int hk = symmetry[k];
                int hl = symmetry[l];
                int kk = k - pitzer_offset[hk];
                int ll = l - pitzer_offset[hl];
                double dik = 0.0;
                double djl = 0.0;
                double dil = 0.0;
                double djk = 0.0;
                if ( hi == hk ) {
                    dik = x->pointer()[d1aoff[hi] + ii * amopi_[hi] + kk];
                }
                if ( hj == hl ) {
                    djl = x->pointer()[d1boff[hj] + jj * amopi_[hj] + ll];
                }
                del2ab->pointer(h)[ij][kl] = x->pointer()[d2aboff[h] + ij * gems_ab[h] + kl] - dik * djl;

            }
        }
    }

    fprintf(rdm_fp,"    v2RDM @||del2(aa)||^2   %20.12lf\n",del2aa->vector_dot(del2aa));
    fprintf(rdm_fp,"    v2RDM @||del2(bb)||^2   %20.12lf\n",del2bb->vector_dot(del2bb));
    fprintf(rdm_fp,"    v2RDM @||del2(ab)||^2   %20.12lf\n",del2ab->vector_dot(del2ab));
    fprintf(rdm_fp,"\n");

    fprintf(rdm_fp,"    v2RDM @Tr[del2(aa)]     %20.12lf\n",del2aa->trace());
    fprintf(rdm_fp,"    v2RDM @Tr[del2(bb)]     %20.12lf\n",del2bb->trace());
    fprintf(rdm_fp,"    v2RDM @Tr[del2(ab)]     %20.12lf\n",del2ab->trace());
    fprintf(rdm_fp,"\n");

    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    ==> v2RDM @D2aa <==\n");
    fprintf(rdm_fp,"\n");

    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,k,l, x_p[d2aaoff[h] + ij * gems_aa[h] + kl]);
                fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",j,i,k,l,-x_p[d2aaoff[h] + ij * gems_aa[h] + kl]);
                fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,l,k,-x_p[d2aaoff[h] + ij * gems_aa[h] + kl]);
                fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",j,i,l,k, x_p[d2aaoff[h] + ij * gems_aa[h] + kl]);
            }
        }
    }
    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    @END\n");

    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    ==> v2RDM @D2bb <==\n");
    fprintf(rdm_fp,"\n");

    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];
            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];
                fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,k,l, x_p[d2bboff[h] + ij * gems_aa[h] + kl]);
                fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",j,j,k,l,-x_p[d2bboff[h] + ij * gems_aa[h] + kl]);
                fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,l,k,-x_p[d2bboff[h] + ij * gems_aa[h] + kl]);
                fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",j,j,l,k, x_p[d2bboff[h] + ij * gems_aa[h] + kl]);
            }
        }
    }
    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    @END\n");

    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    ==> v2RDM @D2ab <==\n");
    fprintf(rdm_fp,"\n");

    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];
            for (int kl = 0; kl < gems_ab[h]; kl++) {
                int k = bas_ab_sym[h][kl][0];
                int l = bas_ab_sym[h][kl][1];
                fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,k,l,x_p[d2aboff[h] + ij * gems_ab[h] + kl]);
            }
        }
    }
    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    @END\n");

    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    ==> v2RDM @D1a <==\n");
    fprintf(rdm_fp,"\n");

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            int ii = i + pitzer_offset[h];
            for (int j = 0; j < amopi_[h]; j++) {
                int jj = j + pitzer_offset[h];
                fprintf(rdm_fp,"%5i %5i %20.12lf\n",i,j,x_p[d1aoff[h] + i * amopi_[h] + j]);
            }
        }
    }
    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    @END\n");

    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    ==> v2RDM @D1b <==\n");
    fprintf(rdm_fp,"\n");

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            int ii = i + pitzer_offset[h];
            for (int j = 0; j < amopi_[h]; j++) {
                int jj = j + pitzer_offset[h];
                fprintf(rdm_fp,"%5i %5i %20.12lf\n",i,j,x_p[d1boff[h] + i * amopi_[h] + j]);
            }
        }
    }
    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    @END\n");

    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    ==> v2RDM @Q1a <==\n");
    fprintf(rdm_fp,"\n");

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            int ii = i + pitzer_offset[h];
            for (int j = 0; j < amopi_[h]; j++) {
                int jj = j + pitzer_offset[h];
                fprintf(rdm_fp,"%5i %5i %20.12lf\n",i,j,x_p[q1aoff[h] + i * amopi_[h] + j]);
            }
        }
    }
    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    @END\n");

    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    ==> v2RDM @Q1b <==\n");
    fprintf(rdm_fp,"\n");

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            int ii = i + pitzer_offset[h];
            for (int j = 0; j < amopi_[h]; j++) {
                int jj = j + pitzer_offset[h];
                fprintf(rdm_fp,"%5i %5i %20.12lf\n",i,j,x_p[q1boff[h] + i * amopi_[h] + j]);
            }
        }
    }
    fprintf(rdm_fp,"\n");
    fprintf(rdm_fp,"    @END\n");

    if ( constrain_q2_ ){

        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    ==> v2RDM @Q2aa <==\n");
        fprintf(rdm_fp,"\n");

        for (int h = 0; h < nirrep_; h++) {
            for (int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                for (int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,k,l, x_p[q2aaoff[h] + ij * gems_aa[h] + kl]);
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",j,i,k,l,-x_p[q2aaoff[h] + ij * gems_aa[h] + kl]);
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,l,k,-x_p[q2aaoff[h] + ij * gems_aa[h] + kl]);
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",j,i,l,k, x_p[q2aaoff[h] + ij * gems_aa[h] + kl]);
                }
            }
        }
        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    @END\n");

        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    ==> v2RDM @Q2bb <==\n");
        fprintf(rdm_fp,"\n");

        for (int h = 0; h < nirrep_; h++) {
            for (int ij = 0; ij < gems_aa[h]; ij++) {
                int i = bas_aa_sym[h][ij][0];
                int j = bas_aa_sym[h][ij][1];
                for (int kl = 0; kl < gems_aa[h]; kl++) {
                    int k = bas_aa_sym[h][kl][0];
                    int l = bas_aa_sym[h][kl][1];
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,k,l, x_p[q2bboff[h] + ij * gems_aa[h] + kl]);
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",j,j,k,l,-x_p[q2bboff[h] + ij * gems_aa[h] + kl]);
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,l,k,-x_p[q2bboff[h] + ij * gems_aa[h] + kl]);
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",j,j,l,k, x_p[q2bboff[h] + ij * gems_aa[h] + kl]);
                }
            }
        }
        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    @END\n");

        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    ==> v2RDM @Q2ab <==\n");
        fprintf(rdm_fp,"\n");

        for (int h = 0; h < nirrep_; h++) {
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,k,l,x_p[q2aboff[h] + ij * gems_ab[h] + kl]);
                }
            }
        }
        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    @END\n");

    }

    if ( constrain_g2_ ){

        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    ==> v2RDM @G2aa/aa <==\n");
        fprintf(rdm_fp,"\n");

        for (int h = 0; h < nirrep_; h++) {
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,k,l,x_p[g2aaoff[h] + (ij             ) * 2 * gems_ab[h] + (kl             )]);
                }
            }
        }
        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    @END\n");

        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    ==> v2RDM @G2aa/bb <==\n");
        fprintf(rdm_fp,"\n");

        for (int h = 0; h < nirrep_; h++) {
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,k,l,x_p[g2aaoff[h] + (ij             ) * 2 * gems_ab[h] + (kl + gems_ab[h])]);
                }
            }
        }
        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    @END\n");

        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    ==> v2RDM @G2bb/aa <==\n");
        fprintf(rdm_fp,"\n");

        for (int h = 0; h < nirrep_; h++) {
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,k,l,x_p[g2aaoff[h] + (ij + gems_ab[h]) * 2 * gems_ab[h] + (kl             )]);
                }
            }
        }
        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    @END\n");

        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    ==> v2RDM @G2bb/bb <==\n");
        fprintf(rdm_fp,"\n");

        for (int h = 0; h < nirrep_; h++) {
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,k,l,x_p[g2aaoff[h] + (ij + gems_ab[h]) * 2 * gems_ab[h] + (kl + gems_ab[h])]);
                }
            }
        }
        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    @END\n");

        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    ==> v2RDM @G2ab <==\n");
        fprintf(rdm_fp,"\n");

        for (int h = 0; h < nirrep_; h++) {
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,k,l,x_p[g2aboff[h] + ij * gems_ab[h] + kl]);
                }
            }
        }
        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    @END\n");
        fprintf(rdm_fp,"\n");

        fprintf(rdm_fp,"    ==> v2RDM @G2ba <==\n");
        fprintf(rdm_fp,"\n");

        for (int h = 0; h < nirrep_; h++) {
            for (int ij = 0; ij < gems_ab[h]; ij++) {
                int i = bas_ab_sym[h][ij][0];
                int j = bas_ab_sym[h][ij][1];
                for (int kl = 0; kl < gems_ab[h]; kl++) {
                    int k = bas_ab_sym[h][kl][0];
                    int l = bas_ab_sym[h][kl][1];
                    fprintf(rdm_fp,"%5i %5i %5i %5i %20.12lf\n",i,j,k,l,x_p[g2baoff[h] + ij * gems_ab[h] + kl]);
                }
            }
        }
        fprintf(rdm_fp,"\n");
        fprintf(rdm_fp,"    @END\n");

    }

    fprintf(rdm_fp,"\n");
}

}
