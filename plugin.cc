/*
 * @BEGIN LICENSE
 *
 * hilbert by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.hpp"

#include "doci/doci_solver.h"
#include "pp2rdm/pp2rdm_solver.h"

using namespace psi;

namespace psi{ namespace hilbert {

extern "C" PSI_API
int read_options(std::string name, Options& options)
{
    if (name == "HILBERT"|| options.read_globals()) {

        /*- SUBSECTION General -*/

        /*- qc solver -*/
        options.add_str("HILBERT_METHOD", "", "DOCI PP2RDM");

        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);

        /*- convergence in the energy -*/
        options.add_double("E_CONVERGENCE", 1e-6);

        /*- convergence in the CI coefficients -*/
        options.add_double("R_CONVERGENCE", 1e-5);

        /*- maximum number of macroiterations -*/
        options.add_int("MAXITER", 50);

        /*- Do write the 2-RDM to disk? All nonzero elements of the 2-RDM will be written.  -*/
        options.add_bool("TPDM_WRITE_FULL",false);

        /*- Do print 2-RDM and 1-RDM to the output file? All nonzero elements of the RDMs will be written.  -*/
        options.add_bool("PRINT_RDMS",false);

        /*- Auxiliary basis set for SCF density fitting computations.
        :ref:`Defaults <apdx:basisFamily>` to a JKFIT basis. -*/
        options.add_str("DF_BASIS_SCF", "");

        /*- What algorithm to use for the SCF computation. See Table :ref:`SCF
        Convergence & Algorithm <table:conv_scf>` for default algorithm for
        different calculation types. -*/
        options.add_str("SCF_TYPE", "DF", "DF CD");

        /*- SUBSECTION DOCI -*/

       /*- maximum size of the subspace in Davidson procedure -*/
        options.add_double("DAVIDSON_MAXDIM", 20);

        /*- Tolerance for Cholesky decomposition of the ERI tensor -*/
        options.add_double("CHOLESKY_TOLERANCE",1e-4);

        /*- Do localize orbitals prior to v2RDM-DOCI? -*/
        options.add_bool("LOCALIZE_ORBITALS",true);

        /*- Do localize virtual orbitals prior to v2RDM-DOCI? -*/
        options.add_bool("LOCALIZE_VIRTUAL_ORBITALS",false);

        /*- Do add random noise to initial orbitals prior to DOCI -*/
        options.add_bool("NOISY_ORBITALS",true);

        /*- Do optimize orbitals? -*/
        options.add_bool("OPTIMIZE_ORBITALS",true);

        /*- SUBSECTION ORBITAL OPTIMIZATION -*/

        /*- algorithm for orbital optimization -*/
        //options.add_str("ORBOPT_ALGORITHM", "HAGER_ZHANG", "STEEPEST_DESCENT HESTENES_STIEFEL DAI_YUAN HAGER_ZHANG KOU_DAI");
        /*- convergence in gradient norm -*/
        options.add_double("ORBOPT_GRADIENT_CONVERGENCE",1.0e-4);
        /*- convergence in energy for rotations -*/
        options.add_double("ORBOPT_ENERGY_CONVERGENCE",1.0e-8);
        /*- flag for using exact expresions for diagonal Hessian element -*/
        options.add_bool("ORBOPT_EXACT_DIAGONAL_HESSIAN",true);
        /*- number of DIIS vectors to keep in orbital optimization -*/
        options.add_int("ORBOPT_NUM_DIIS_VECTORS",0);
        /*- maximum number of iterations for orbital optimization -*/
        options.add_int("ORBOPT_MAXITER",1);
        /*- Do write a ORBOPT output file?  If so, the filename will end in
        .molden, and the prefix is determined by |globals__writer_file_label|
        (if set), or else by the name of the output file plus the name of
        the current molecule. -*/
        options.add_bool("ORBOPT_WRITE", false);

        /*- SUBSECTION pp2RDM -*/

        /*- do check analytic gradient for accuracy? -*/
        options.add_bool("CHECK_GRADIENT",false);

        /*- do check analytic hessian for accuracy? -*/
        options.add_bool("CHECK_HESSIAN",false);

        /*- restart from checkpoint file? -*/
        options.add_bool("RESTART_FROM_CHECKPOINT_FILE",false);

        /*- algorithm type -*/
        options.add_str("P2RDM_ALGORITHM","PROJECTION","PROJECTION LBFGS NEWTON_RAPHSON");

        /*- Which parametric 2-RDM method is called? -*/
        options.add_str("P2RDM_TYPE","K","K CEPA(0) CEPA(1) CID ACPF AQCC CCD");

        /*- Do print 1- and 2-electron to the output file? Only J-, K-, and L-type integrals will be printed. -*/
        options.add_bool("PRINT_INTEGRALS",false);

    }

    return true;
}

extern "C" PSI_API
SharedWavefunction hilbert(SharedWavefunction ref_wfn, Options& options)
{

    if ( options.get_str("HILBERT_METHOD") == "DOCI") {

        std::shared_ptr<doci::DOCISolver> doci (new doci::DOCISolver(ref_wfn,options));
        double energy = doci->compute_energy();
        return (std::shared_ptr<Wavefunction>)doci;

    }else if ( options.get_str("HILBERT_METHOD") == "PP2RDM") {

        std::shared_ptr<pp2rdm::pp2RDMSolver> pp2rdm (new pp2rdm::pp2RDMSolver(ref_wfn,options));
        double energy = pp2rdm->compute_energy();
        return (std::shared_ptr<Wavefunction>)pp2rdm;

    }
}

}} // End namespaces

