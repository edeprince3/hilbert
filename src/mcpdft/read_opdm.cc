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

#include "psi4/psi4-dec.h"
#include <psi4/psifiles.h>
#include <psi4/libiwl/iwl.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libmints/matrix.h>
#include "psi4/libmints/vector.h"
#include <psi4/libpsi4util/PsiOutStream.h>

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libqt/qt.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "mcpdft_solver.h"

using namespace psi;
using namespace std;

namespace hilbert{

void MCPDFTSolver::ReadOPDM() {

    std::shared_ptr<PSIO> psio (new PSIO());

    // TODO: should be added back when reading-in the density PSIOH files 
    // psio->set_pid("18332");

    if ( !psio->exists(PSIF_V2RDM_D1A) ) throw PsiException("No D1a on disk",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D1B) ) throw PsiException("No D1b on disk",__FILE__,__LINE__);

    // D1a

    psio->open(PSIF_V2RDM_D1A,PSIO_OPEN_OLD);

    long int na;
    psio->read_entry(PSIF_V2RDM_D1A,"length",(char*)&na,sizeof(long int));

    opdm_a_ = (opdm *)malloc(na * sizeof(opdm));
    psio->read_entry(PSIF_V2RDM_D1A,"D1a",(char*)opdm_a_,na * sizeof(opdm));
    psio->close(PSIF_V2RDM_D1A,1);

    for (int n = 0; n < na; n++) {

        int i = opdm_a_[n].i;
        int j = opdm_a_[n].j;

        int hi = symmetry_[i];
        int hj = symmetry_[j];

        if ( hi != hj ) {
            throw PsiException("error: something is wrong with the symmetry of the alpha OPDM",__FILE__,__LINE__);
        }

        int ii = i - pitzer_offset_[hi];
        int jj = j - pitzer_offset_[hi];

        Da_->pointer(hi)[ii][jj] = opdm_a_[n].value;

    }

    // D1b

    psio->open(PSIF_V2RDM_D1B,PSIO_OPEN_OLD);

    long int nb;
    psio->read_entry(PSIF_V2RDM_D1B,"length",(char*)&nb,sizeof(long int));

    opdm_b_ = (opdm *)malloc(nb * sizeof(opdm));
    psio->read_entry(PSIF_V2RDM_D1B,"D1b",(char*)opdm_b_,nb * sizeof(opdm));
    psio->close(PSIF_V2RDM_D1B,1);

   for (int n = 0; n < nb; n++) {

        int i = opdm_b_[n].i;
        int j = opdm_b_[n].j;

        int hi = symmetry_[i];
        int hj = symmetry_[j];

        if ( hi != hj ) {
            throw PsiException("error: something is wrong with the symmetry of the beta OPDM",__FILE__,__LINE__);
        }

        int ii = i - pitzer_offset_[hi];
        int jj = j - pitzer_offset_[hi];

        Db_->pointer(hi)[ii][jj] = opdm_b_[n].value;

    }

    if ( !is_low_memory_ ) {
        outfile->Printf("\n");
        outfile->Printf("    ==> Build Rho's ...\n");
        BuildRhoFast(na,nb);
        outfile->Printf("    ... Done. <==\n\n");
    }
    free(opdm_a_);
    free(opdm_b_);
}

void MCPDFTSolver::ReadCIOPDM(std::shared_ptr<Matrix> D, const char* fileName) {

    std::ifstream dataIn;

    dataIn.open(fileName);

    if (!dataIn) throw PsiException("No D1 on disk",__FILE__,__LINE__);
    else {
        double ** dp = D->pointer();
        for (int i = 0; i < nmo_; i++) {
            for (int j = 0; j < nmo_; j++) {
                dataIn >> dp[i][j];
                if (dp[i][j] < 1.0e-20)
                    dp[i][j] = 0.0;
            }
        }
        dataIn.close();
    }
}

} //end namespaces


