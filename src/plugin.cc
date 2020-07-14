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
#include <psi4/libpsi4util/PsiOutStream.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libpsio/psio.hpp>

#include <v2rdm_doci/v2rdm_solver.h>
#include <doci/doci_solver.h>
#include <pp2rdm/pp2rdm_solver.h>

using namespace psi;

namespace hilbert {

extern "C" PSI_API
int read_options(std::string name, Options& options)
{
    if (name == "HILBERT"|| options.read_globals()) {

        /*- SUBSECTION General -*/

        /*- qc solver -*/
        options.add_str("HILBERT_METHOD", "", "DOCI PP2RDM V2RDM_DOCI");

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

        /*- algorithm for orbital optimization. only valid for v2rdm-doci -*/
        options.add_str("ORBOPT_ALGORITHM", "HAGER_ZHANG", "STEEPEST_DESCENT HESTENES_STIEFEL DAI_YUAN HAGER_ZHANG KOU_DAI");

        /*- frequency of orbital optimization.  optimization occurs every 
        orbopt_frequency iterations. only valid for v2rdm-doci -*/
        options.add_int("ORBOPT_FREQUENCY",500);

        /*- convergence in gradient norm -*/
        options.add_double("ORBOPT_GRADIENT_CONVERGENCE",1.0e-4);

        /*- convergence in energy for rotations -*/
        options.add_double("ORBOPT_ENERGY_CONVERGENCE",1.0e-8);

        /*- do rotate active-active orbital pairs. No methods in Hilbert use this flag currently !expert -*/
        options.add_bool("ORBOPT_ACTIVE_ACTIVE_ROTATIONS",false);

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

        /*- Which parametric 2-RDM method is called? Set by driver. !expert -*/
        options.add_str("P2RDM_TYPE","K","K CEPA(0) CEPA(1) CID ACPF AQCC CCD");

        /*- Do print 1- and 2-electron to the output file? Only J-, K-, and L-type integrals will be printed. -*/
        options.add_bool("PRINT_INTEGRALS",false);

        /*- SUBSECTION v2RDM-DOCI -*/

        /* Do v2RDM-DOCI gradient? !expert */
        options.add_str("DERTYPE", "NONE", "NONE FIRST");

        /*- Do semicanonicalize orbitals? -*/
        options.add_bool("SEMICANONICALIZE_ORBITALS",false);

        /*- Type of guess -*/
        options.add_str("TPDM_GUESS","RANDOM", "RANDOM HF");

        /*- Do save progress in a checkpoint file? -*/
        options.add_bool("WRITE_CHECKPOINT_FILE",false);

        /*- Frequency of checkpoint file generation.  The checkpoint file is 
        updated every CHECKPOINT_FREQUENCY iterations.  The default frequency
        will be ORBOPT_FREQUENCY. -*/
        options.add_int("CHECKPOINT_FREQUENCY",500);

        /*- Frequency with which the pentalty-parameter, mu, is updated. mu is
        updated every MU_UPDATE_FREQUENCY iterations.   -*/
        options.add_int("MU_UPDATE_FREQUENCY",1000);

        /*- The type of 2-positivity computation -*/
        options.add_str("POSITIVITY", "DQG", "DQG D DQ DG DQGT2 DQGT1 DQGT1T2 DQGT2");

        /*- Do constrain D3/D2 mapping? -*/
        options.add_bool("CONSTRAIN_D3",false);

        /*- convergence for conjugate gradient solver. currently not used. -*/
        options.add_double("CG_CONVERGENCE", 1e-9);

        /*- maximum number of conjugate gradient iterations -*/
        options.add_int("CG_MAXITER", 10000);

    }

    return true;
}

extern "C" PSI_API
SharedWavefunction hilbert(SharedWavefunction ref_wfn, Options& options)
{

    if ( options.get_str("HILBERT_METHOD") == "DOCI") {

        std::shared_ptr<DOCISolver> doci (new DOCISolver(ref_wfn,options));
        double energy = doci->compute_energy();
        return (std::shared_ptr<Wavefunction>)doci;

    }else if ( options.get_str("HILBERT_METHOD") == "PP2RDM") {

        std::shared_ptr<pp2RDMSolver> pp2rdm (new pp2RDMSolver(ref_wfn,options));
        double energy = pp2rdm->compute_energy();
        return (std::shared_ptr<Wavefunction>)pp2rdm;

    }else if ( options.get_str("HILBERT_METHOD") == "V2RDM_DOCI") {

        std::shared_ptr<v2RDMSolver> v2rdm_doci (new v2RDMSolver(ref_wfn,options));
        double energy = v2rdm_doci->compute_energy();
        return (std::shared_ptr<Wavefunction>)v2rdm_doci;

    }

    return ref_wfn;
}

} // End namespaces

