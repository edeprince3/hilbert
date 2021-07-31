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

#include <psi4/libmints/mintshelper.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libtrans/mospace.h>

#include "v2rdm_solver.h"

#include <misc/omp.h>
#include <misc/threeindexintegrals.h>

using namespace psi;

namespace hilbert{

void v2RDMSolver::ComputeNaturalOrbitals() {

    //if ( nalpha_ != nbeta_ || !constrain_spin_ ) {
    //    throw PsiException("natural orbital transformation only implemented for closed-shell systems with spin symmetry enforced",__FILE__,__LINE__);
    //}

    if ( nalpha_ != nbeta_ || !constrain_spin_ ) {

        outfile->Printf("\n");
        outfile->Printf("    ==> Warning <==\n");
        outfile->Printf("\n");
        outfile->Printf("        The natural orbitals computed here for EKT and FCIDUMP are\n");
        outfile->Printf("        the natural orbitals of the spin-free 1-RDM. On the other\n");
        outfile->Printf("        hand, the natural orbital occupation numbers printed below\n");
        outfile->Printf("        correspond to natural spin orbitals.\n");
        outfile->Printf("\n");

    }

    SharedMatrix D (new Matrix(nirrep_,amopi_,amopi_));
    SharedMatrix eigvec (new Matrix(nirrep_,amopi_,amopi_));
    SharedVector eigval (new Vector("Natural Orbital Occupation Numbers",nirrep_,amopi_));

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                D->pointer(h)[i][j]  = x->pointer()[d1aoff[h]+i*amopi_[h]+j];
                D->pointer(h)[i][j] += x->pointer()[d1boff[h]+i*amopi_[h]+j];
            }
        }
    }
    D->diagonalize(eigvec,eigval,descending);

    // build AO/NO transformation matrix (both Ca_ and Cb_)
    for (int h = 0; h < nirrep_; h++) {

        for (int mu = 0; mu < nsopi_[h]; mu++) {

            double *  temp = (double*)malloc(nmopi_[h]*sizeof(double));
            double ** ep   = eigvec->pointer(h);

            // Ca_

            double ** cp = Ca_->pointer(h);
            for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
                double dum = 0.0;
                for (int j = rstcpi_[h] + frzcpi_[h]; j < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; j++) {
                    dum += cp[mu][j] * ep[j-rstcpi_[h]-frzcpi_[h]][i-rstcpi_[h]-frzcpi_[h]];
                }
                temp[i] = dum;
            }
            for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
                cp[mu][i] = temp[i];
            }

/*
            // Cb_
            cp = Cb_->pointer(h);
            for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
                double dum = 0.0;
                for (int j = rstcpi_[h] + frzcpi_[h]; j < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; j++) {
                    dum += cp[mu][j] * ep[j-rstcpi_[h]-frzcpi_[h]][i-rstcpi_[h]-frzcpi_[h]];
                }
                temp[i] = dum;
            }
            for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
                cp[mu][i] = temp[i];
            }
*/

            free(temp);
        }
    }

    // transform the alpha 1-RDM to natural orbital basis
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                D->pointer(h)[i][j]  = x->pointer()[d1aoff[h]+i*amopi_[h]+j];
            }
        }
    }
    D->transform(eigvec);
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                x->pointer()[d1aoff[h]+i*amopi_[h]+j] = D->pointer(h)[i][j];
            }
        }
    }

    // transform the beta 1-RDM to natural orbital basis
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                D->pointer(h)[i][j]  = x->pointer()[d1boff[h]+i*amopi_[h]+j];
            }
        }
    }
    D->transform(eigvec);
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < amopi_[h]; i++) {
            for (int j = 0; j < amopi_[h]; j++) {
                x->pointer()[d1boff[h]+i*amopi_[h]+j] = D->pointer(h)[i][j];
            }
        }
    }


    // transform the ab block of the 2-RDM to natural orbital basis
    std::shared_ptr<Vector> tempx (new Vector(n_primal_));

    // D2ab(ij,kl): transform index 4
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];
    
                            double dum = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                int kp = ibas_ab_sym[h][k][p];

                                dum += x->pointer()[d2aboff[h]+ij*gems_ab[h]+kp] * ep[pp][ll];

                            }

                            int kl = ibas_ab_sym[h][k][l];

                            tempx->pointer()[d2aboff[h]+ij*gems_ab[h]+kl] = dum;

                        }
                    }
                }
            }
        }
    }

    // D2ab(ij,lk): transform index 3
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];
    
                            double dum = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                int pk = ibas_ab_sym[h][p][k];

                                dum += tempx->pointer()[d2aboff[h]+ij*gems_ab[h]+pk] * ep[pp][ll];

                            }

                            int lk = ibas_ab_sym[h][l][k];

                            x->pointer()[d2aboff[h]+ij*gems_ab[h]+lk] = dum;

                        }
                    }
                }
            }
        }
    }

    // D2ab(kl,ij): transform index 2
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];
    
                            double dum = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                int kp = ibas_ab_sym[h][k][p];

                                dum += x->pointer()[d2aboff[h]+kp*gems_ab[h]+ij] * ep[pp][ll];

                            }

                            int kl = ibas_ab_sym[h][k][l];

                            tempx->pointer()[d2aboff[h]+kl*gems_ab[h]+ij] = dum;

                        }
                    }
                }
            }
        }
    }

    // D2ab(lk,ij): transform index 1
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];
    
                            double dum = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                int pk = ibas_ab_sym[h][p][k];

                                dum += tempx->pointer()[d2aboff[h]+pk*gems_ab[h]+ij] * ep[pp][ll];

                            }

                            int lk = ibas_ab_sym[h][l][k];

                            x->pointer()[d2aboff[h]+lk*gems_ab[h]+ij] = dum;

                        }
                    }
                }
            }
        }
    }

    // now transform the aa and bb blocks of the 2-RDM to natural orbital basis

    // first, unpack D2aa/bb into full antisymmetrized D2aa/D2bb
    std::shared_ptr<Vector> tempx2 (new Vector(n_primal_));
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];

            int ij_ab = ibas_ab_sym[h][i][j];
            int ji_ab = ibas_ab_sym[h][j][i];

            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];

                int kl_ab = ibas_ab_sym[h][k][l];
                int lk_ab = ibas_ab_sym[h][l][k];

                tempx2->pointer()[d2aboff[h]+ij_ab*gems_ab[h]+kl_ab] =  x->pointer()[d2aaoff[h] + ij*gems_aa[h] + kl];
                tempx2->pointer()[d2aboff[h]+ji_ab*gems_ab[h]+kl_ab] = -x->pointer()[d2aaoff[h] + ij*gems_aa[h] + kl];
                tempx2->pointer()[d2aboff[h]+ij_ab*gems_ab[h]+lk_ab] = -x->pointer()[d2aaoff[h] + ij*gems_aa[h] + kl];
                tempx2->pointer()[d2aboff[h]+ji_ab*gems_ab[h]+lk_ab] =  x->pointer()[d2aaoff[h] + ij*gems_aa[h] + kl];

                tempx2->pointer()[q2aboff[h]+ij_ab*gems_ab[h]+kl_ab] =  x->pointer()[d2bboff[h] + ij*gems_aa[h] + kl];
                tempx2->pointer()[q2aboff[h]+ji_ab*gems_ab[h]+kl_ab] = -x->pointer()[d2bboff[h] + ij*gems_aa[h] + kl];
                tempx2->pointer()[q2aboff[h]+ij_ab*gems_ab[h]+lk_ab] = -x->pointer()[d2bboff[h] + ij*gems_aa[h] + kl];
                tempx2->pointer()[q2aboff[h]+ji_ab*gems_ab[h]+lk_ab] =  x->pointer()[d2bboff[h] + ij*gems_aa[h] + kl];
            }
        }
    }

    // D2aa(ij,kl): transform index 4
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];
    
                            double duma = 0.0;
                            double dumb = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                int kp = ibas_ab_sym[h][k][p];

                                duma += tempx2->pointer()[d2aboff[h]+ij*gems_ab[h]+kp] * ep[pp][ll];
                                dumb += tempx2->pointer()[q2aboff[h]+ij*gems_ab[h]+kp] * ep[pp][ll];

                            }

                            int kl = ibas_ab_sym[h][k][l];

                            tempx->pointer()[d2aboff[h]+ij*gems_ab[h]+kl] = duma;
                            tempx->pointer()[q2aboff[h]+ij*gems_ab[h]+kl] = dumb;

                        }
                    }
                }
            }
        }
    }

    // D2aa(ij,lk): transform index 3
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];
    
                            double duma = 0.0;
                            double dumb = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                int pk = ibas_ab_sym[h][p][k];

                                duma += tempx->pointer()[d2aboff[h]+ij*gems_ab[h]+pk] * ep[pp][ll];
                                dumb += tempx->pointer()[q2aboff[h]+ij*gems_ab[h]+pk] * ep[pp][ll];

                            }

                            int lk = ibas_ab_sym[h][l][k];

                            tempx2->pointer()[d2aboff[h]+ij*gems_ab[h]+lk] = duma;
                            tempx2->pointer()[q2aboff[h]+ij*gems_ab[h]+lk] = dumb;

                        }
                    }
                }
            }
        }
    }

    // D2aa(kl,ij): transform index 2
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];
    
                            double duma = 0.0;
                            double dumb = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                int kp = ibas_ab_sym[h][k][p];

                                duma += tempx2->pointer()[d2aboff[h]+kp*gems_ab[h]+ij] * ep[pp][ll];
                                dumb += tempx2->pointer()[q2aboff[h]+kp*gems_ab[h]+ij] * ep[pp][ll];

                            }

                            int kl = ibas_ab_sym[h][k][l];

                            tempx->pointer()[d2aboff[h]+kl*gems_ab[h]+ij] = duma;
                            tempx->pointer()[q2aboff[h]+kl*gems_ab[h]+ij] = dumb;

                        }
                    }
                }
            }
        }
    }

    // D2aa(lk,ij): transform index 1
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_ab[h]; ij++) {

            for (int hk = 0; hk < nirrep_; hk++) {
                for (int kk = 0; kk < amopi_[hk]; kk++) {
                    int k = kk + pitzer_offset[hk];

                    for (int hl = 0; hl < nirrep_; hl++) {

                        int hkl = SymmetryPair(hk,hl);
                        if ( hkl != h ) continue;

                        double ** ep   = eigvec->pointer(hl);

                        for (int ll = 0; ll < amopi_[hl]; ll++) {
                            int l = ll + pitzer_offset[hl];
    
                            double duma = 0.0;
                            double dumb = 0.0;
                            for (int pp = 0; pp < amopi_[hl]; pp++) {

                                int p = pp + pitzer_offset[hl];

                                int pk = ibas_ab_sym[h][p][k];

                                duma += tempx->pointer()[d2aboff[h]+pk*gems_ab[h]+ij] * ep[pp][ll];
                                dumb += tempx->pointer()[q2aboff[h]+pk*gems_ab[h]+ij] * ep[pp][ll];

                            }

                            int lk = ibas_ab_sym[h][l][k];

                            tempx2->pointer()[d2aboff[h]+lk*gems_ab[h]+ij] = duma;
                            tempx2->pointer()[q2aboff[h]+lk*gems_ab[h]+ij] = dumb;

                        }
                    }
                }
            }
        }
    }

    // lastly, repack D2aa/bb from full antisymmetrized D2aa/D2bb
    for (int h = 0; h < nirrep_; h++) {
        for (int ij = 0; ij < gems_aa[h]; ij++) {
            int i = bas_aa_sym[h][ij][0];
            int j = bas_aa_sym[h][ij][1];

            int ij_ab = ibas_ab_sym[h][i][j];
            int ji_ab = ibas_ab_sym[h][j][i];

            for (int kl = 0; kl < gems_aa[h]; kl++) {
                int k = bas_aa_sym[h][kl][0];
                int l = bas_aa_sym[h][kl][1];

                int kl_ab = ibas_ab_sym[h][k][l];
                int lk_ab = ibas_ab_sym[h][l][k];

                x->pointer()[d2aaoff[h] + ij*gems_aa[h] + kl] = tempx2->pointer()[d2aboff[h]+ij_ab*gems_ab[h]+kl_ab];
                x->pointer()[d2bboff[h] + ij*gems_aa[h] + kl] = tempx2->pointer()[q2aboff[h]+ij_ab*gems_ab[h]+kl_ab];

            }
        }
    }

    // if dumping integrals to disk, transform these to the NO basis as well.

    if ( options_.get_bool("FCIDUMP") || options_.get_bool("EXTENDED_KOOPMANS") ) {

        if ( is_df_ ) {
            free(Qmo_);

            double start = omp_get_wtime();
            ThreeIndexIntegrals(reference_wavefunction_,nQ_,memory_);

            Qmo_ = (double*)malloc(nmo_*(nmo_+1)/2*nQ_*sizeof(double));
            memset((void*)Qmo_,'\0',nmo_*(nmo_+1)/2*nQ_*sizeof(double));

            std::shared_ptr<PSIO> psio(new PSIO());
            psio->open(PSIF_DCC_QMO,PSIO_OPEN_OLD);
            psio->read_entry(PSIF_DCC_QMO,"(Q|mn) Integrals",(char*)Qmo_,sizeof(double)*nQ_ * nmo_*(nmo_+1)/2);
            psio->close(PSIF_DCC_QMO,1);

        }else {
            std::vector<std::shared_ptr<MOSpace> > spaces;
            spaces.push_back(MOSpace::all);
            std::shared_ptr<IntegralTransform> ints(new IntegralTransform(reference_wavefunction_, spaces,
                IntegralTransform::TransformationType::Restricted, IntegralTransform::OutputType::IWLOnly,
                IntegralTransform::MOOrdering::PitzerOrder, IntegralTransform::FrozenOrbitals::None, false));
            ints->set_dpd_id(0);
            ints->set_keep_iwl_so_ints(true);
            ints->set_keep_dpd_so_ints(true);
            ints->initialize();
            ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
            GetTEIFromDisk();
        }
    }

}

void v2RDMSolver::PrintNaturalOrbitalOccupations() {

    SharedMatrix Da (new Matrix(nirrep_,nmopi_,nmopi_));
    SharedMatrix eigveca (new Matrix(nirrep_,nmopi_,nmopi_));
    SharedVector eigvala (new Vector("Natural Orbital Occupation Numbers (alpha)",nirrep_,nmopi_));

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h] + rstcpi_[h]; i++) {
            Da->pointer(h)[i][i] = 1.0;
        }
        for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h] - rstvpi_[h] - frzvpi_[h]; i++) {
            for (int j = rstcpi_[h] + frzcpi_[h]; j < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; j++) {
                Da->pointer(h)[i][j] = x->pointer()[d1aoff[h]+(i-rstcpi_[h]-frzcpi_[h])*amopi_[h]+(j-rstcpi_[h]-frzcpi_[h])];
            }
        }
    }
    SharedMatrix saveda ( new Matrix(Da) );
    Da->diagonalize(eigveca,eigvala,descending);
    eigvala->print();

    SharedMatrix Db (new Matrix(nirrep_,nmopi_,nmopi_));
    SharedMatrix eigvecb (new Matrix(nirrep_,nmopi_,nmopi_));
    SharedVector eigvalb (new Vector("Natural Orbital Occupation Numbers (beta)",nirrep_,nmopi_));
    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < rstcpi_[h] + frzcpi_[h]; i++) {
            Db->pointer(h)[i][i] = 1.0;
        }
        for (int i = rstcpi_[h] + frzcpi_[h]; i < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; i++) {
            for (int j = rstcpi_[h] + frzcpi_[h]; j < nmopi_[h]-rstvpi_[h]-frzvpi_[h]; j++) {
                Db->pointer(h)[i][j] = x->pointer()[d1boff[h]+(i-rstcpi_[h]-frzcpi_[h])*amopi_[h]+(j-rstcpi_[h]-frzcpi_[h])];
            }
        }
    }
    Db->diagonalize(eigvecb,eigvalb,descending);
    eigvalb->print();

}


} // End namespaces

