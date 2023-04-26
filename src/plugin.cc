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
#include <psi4/libpsi4util/process.h>

#include <v2rdm_casscf/v2rdm_solver.h>
#include <v2rdm_doci/v2rdm_doci_solver.h>
#include <doci/doci_solver.h>
#include <pp2rdm/pp2rdm_solver.h>
#include <p2rdm/p2rdm_solver.h>
#include <jellium/jellium_scf_solver.h>
#include <polaritonic_scf/rhf.h>
#include <polaritonic_scf/uhf.h>
#include <polaritonic_scf/rohf.h>
#include <polaritonic_scf/rks.h>
#include <polaritonic_scf/uks.h>
#include <polaritonic_scf/uccsd.h>
#include <polaritonic_scf/rcis.h>
#include <polaritonic_scf/rtddft.h>
#include <polaritonic_scf/utddft.h>
#include <misc/backtransform_tpdm.h>

#ifdef USE_QED_CC
    #include "python_api/python_helpers.h"
    #include <tiledarray.h>
    #include "cc_cavity/include/cc_cavity.h"
    #include "cc_cavity/include/derived/qed_ccsd.h"
    #if MAX_CC_LEVEL >= 3
        #include "cc_cavity/include/derived/qed_ccsdt.h"
    #endif
    #if MAX_CC_LEVEL >= 4
        #include "cc_cavity/include/derived/qed_ccsdtq.h"
    #endif
#endif

#ifdef USE_QED_EOM_CC
    #include "cc_cavity/include/derived/eom_ee_ccsd.h"
    #include "cc_cavity/include/derived/eom_ee_qed_ccsd.h"
    #include "cc_cavity/include/derived/eom_ea_driver.h"
    #include "cc_cavity/include/derived/eom_ee_rdm.h"
#endif


using namespace psi;

namespace hilbert {

extern "C" PSI_API
int read_options(std::string name, Options& options)
{
    if (name == "HILBERT"|| options.read_globals()) {

        /*- SUBSECTION General -*/

        /*- qc solver. used internally !expert -*/
        options.add_str("HILBERT_METHOD", "", "DOCI P2RDM PP2RDM V2RDM_DOCI V2RDM_CASSCF JELLIUM_SCF POLARITONIC_RHF POLARITONIC_UHF POLARITONIC_ROHF POLARITONIC_UKS POLARITONIC_RKS POLARITONIC_RCIS POLARITONIC_UCCSD POLARITONIC_TDDFT CC_CAVITY");


        /*- Do DIIS? -*/
        options.add_bool("DIIS", true);

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

        /*- Do print t-and z-amplitudes to the output file? -*/
        options.add_bool("PRINT_PCCD_AMPLITUDES",false);

        /*- Auxiliary basis set for SCF density fitting computations.
        :ref:`Defaults <apdx:basisFamily>` to a JKFIT basis. -*/
        options.add_str("DF_BASIS_SCF", "");

        /*- What algorithm to use for the SCF computation. See Table :ref:`SCF
        Convergence & Algorithm <table:conv_scf>` for default algorithm for
        different calculation types. -*/
        options.add_str("SCF_TYPE", "DF", "DF CD");

        /*- SUBSECTION DOCI -*/

       /*- maximum size of Davidson subspace (will be multiplied by number of desired roots) -*/
        options.add_double("DAVIDSON_MAXDIM", 20);

        /*- Tolerance for Cholesky decomposition of the ERI tensor -*/
        options.add_double("CHOLESKY_TOLERANCE",1e-4);

        /*- Do localize orbitals prior to v2RDM-DOCI? -*/
        options.add_bool("LOCALIZE_ORBITALS",false);

        /*- Do localize virtual orbitals prior to v2RDM-DOCI? -*/
        options.add_bool("LOCALIZE_VIRTUAL_ORBITALS",false);

        /*- Do add random noise to initial orbitals prior to DOCI -*/
        options.add_bool("NOISY_ORBITALS",false);

        /*- Do optimize orbitals? -*/
        options.add_bool("OPTIMIZE_ORBITALS",true);

        /*- SUBSECTION ORBITAL OPTIMIZATION -*/

        options.add_bool("MOLDEN_WRITE", false);
        /*- Do write a MOLDEN file for guess orbitals?  If so, the filename will
        end in .guess.molden, and the prefix is determined by 
        |globals__writer_file_label| (if set), or else by the name of the output
        file plus the name of the current molecule. -*/

        options.add_bool("GUESS_ORBITALS_WRITE", false);
        /*- Do write a ORBOPT output file?  If so, the filename will end in
        .molden, and the prefix is determined by |globals__writer_file_label|
        (if set), or else by the name of the output file plus the name of
        the current molecule. -*/

        /*- flag to optimize orbitals using a one-step type approach -*/
        options.add_bool("ORBOPT_ONE_STEP",true);

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
        options.add_int("ORBOPT_MAXITER",10);

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

        /*- File containing previous primal/dual solutions and integrals. -*/
        options.add_str("RESTART_FROM_CHECKPOINT_FILE","");

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
        options.add_str("POSITIVITY", "DQG", "DQG D DQ DG DQGT1 DQGT2 DQGT1T2 3POS");

        /*- Do enforce generalized pauli constraints -*/ 
        options.add_str("GPC_CONSTRAINTS","NONE", "NONE 1RDM 2RDM");

        /*- Do constrain D3 to D2 mapping? -*/
        options.add_bool("CONSTRAIN_D3",false);

        /*- Do constrain E3 to D3 mapping? -*/
        options.add_bool("CONSTRAIN_E3",false);

        /*- Do constrain F3 to D3 mapping? -*/
        options.add_bool("CONSTRAIN_F3",false);

        /*- Do constrain Q3 to D3 mapping? -*/
        options.add_bool("CONSTRAIN_Q3",false);

        /*- convergence for conjugate gradient solver. currently not used. -*/
        options.add_double("CG_CONVERGENCE", 1e-9);

        /*- maximum number of conjugate gradient iterations -*/
        options.add_int("CG_MAXITER", 10000);

        /*- SUBSECTION v2RDM-CASSCF -*/

        /*- SDP solver -*/
        options.add_str("SDP_SOLVER","BPSDP", "BPSDP RRSDP");

        /*- do use hubbard model? -*/
        options.add_bool("HUBBARD_HAMILTONIAN",false);

        /*- hubbard hopping integral -*/
        options.add_double("HUBBARD_T",1.0);

        /*- hubbard on-site repulsion -*/
        options.add_double("HUBBARD_U",1.0);

        /*- number of sites in hubbard model -*/
        options.add_int("N_HUBBARD_SITES",4);

        /*- total number of spins in hubbard model -*/
        options.add_int("N_HUBBARD_SPINS",4);

        /*- multiplicity in hubbard model -*/
        options.add_int("HUBBARD_MULTIPLICITY",1);

        /*- fractional charge -*/
        options.add_double("FRACTIONAL_CHARGE", 0.0);

        /*- do extended koopmans theorem computation? -*/
        options.add_bool("EXTENDED_KOOPMANS",false);

        /*- Do v2RDM-CASSCF gradient? !expert -*/
        options.add_str("DERTYPE", "NONE", "NONE FIRST");

        /* Do write fcidump files? -*/
        options.add_bool("FCIDUMP", false);

        /*- Rotate guess orbitals -*/
        options.add("MCSCF_ROTATE", new ArrayType());

        /*- Do compute natural orbitals and transform 1- and 2-RDM to the natural orbital basis? 
        The OPDM and Ca/Cb matrices pushed onto the wavefunction will correspond to the natural orbital basis -*/
        options.add_bool("NAT_ORBS",false);

        /*- Do write the 1-RDM to disk? All nonzero elements of the 1-RDM will be written.  -*/
        options.add_bool("OPDM_WRITE_FULL",false);

        /*- Do write the spin-free 2-RDM to disk? All nonzero elements of the 2-RDM will be written.  -*/
        options.add_bool("TPDM_WRITE_SPIN_FREE",false);

        /*- Do write the 2-RDM to disk? Only the nonzero elements of the active 2-RDM will be written. -*/
        options.add_bool("TPDM_WRITE",false);

        /*- Do write the 3-RDM to disk? -*/
        options.add_bool("3PDM_WRITE",false);

        /*- A parameter introduced by Mazziotti [PRL 106, 083001 (2011)] to "increase the
        sensitivity of y on the deviation of x from primal feasibility."  Should 
        lie on the interval [1.0, 1.6]. -*/
        options.add_double("TAU_PARAMETER",1.0);

        /*- Do constrain D4 to D3 mapping? -*/
        options.add_bool("CONSTRAIN_D4",false);

        /*- Do constrain spin squared? -*/
        options.add_bool("CONSTRAIN_SPIN", true);

        /*- Do constrain sz? -*/
        options.add_bool("CONSTRAIN_SZ", true);

        /*- SUBSECTION JELLIUM -*/

        /*- An array containing the number of doubly-occupied orbitals per irrep
        (in Cotton order) -*/
        options.add("DOCC", new ArrayType());

        /*- The length of the box in nm. No default. If not specified, the
        box length is chosen to satisfy <rho> = 1e-/a0^3 -*/
        options.add_double("JELLIUM_BOX_LENGTH",1.0);

        /*- The number of grid points for the Gauss-Legendre quadrature -*/
        options.add_int("N_GRID_POINTS", 10);

        /*- The number of electrons -*/
        options.add_int("N_ELECTRONS", 2);

        /*- The number of basis functions -*/
        options.add_int("N_BASIS_FUNCTIONS", 26);

        /*- The length of the box in nm -*/
        options.add_double("LENGTH", 1.0);
        //options.add_double("LENGTH", 0.166245);

        ///*- The density of the box in e/nm^3 -*/
        //options.add_double("DENSITY", 92);    

        /*- The number of electronic states to computed, per irreducible
        representation -*/
        options.add("ROOTS_PER_IRREP", new ArrayType());

        /*- Do smart guess in Davidson? Requires exact hamiltonian elements 
        and could get expensive, default = false -*/
        options.add_bool("JELLIUM_CIS_SMART_GUESS", false);

        /*- SUBSECTION POLARITONIC SCF -*/

        /*- functional for cavity QED-DFT -*/
        options.add_str("CAVITY_QED_DFT_FUNCTIONAL", "B3LYP");

        /*- number of photon number states -*/
        options.add_int("N_PHOTON_STATES", 1);

        /*- do use coherent-state basis? !expert -*/
        options.add_bool("USE_COHERENT_STATE_BASIS", true);

        /*- cavity excitation energy for the modes along the x, y and z axis (a.u.) -*/
        options.add("CAVITY_FREQUENCY",new ArrayType());

        /*- cavity coupling strength (a.u.) -*/
        options.add("CAVITY_COUPLING_STRENGTH",new ArrayType());

        /*- do include u0 in polaritioinic ccsd? -*/
        options.add_bool("POLARITONIC_CC_INCLUDE_U0",false);

        /*- do include u1 in polaritioinic ccsd? -*/
        options.add_bool("POLARITONIC_CC_INCLUDE_U1",false);

        /*- do include u2 in polaritioinic ccsd? -*/
        options.add_bool("POLARITONIC_CC_INCLUDE_U2",false);

        /*- do use TDA in TDDFT? -*/
        options.add_bool("TDSCF_TDA",false);

        /*- do relax orbitals in QED-SCF [unlike QED-TDDFT described in J. Chem. Phys. 155, 064107 (2021)?] -*/
        options.add_bool("QED_USE_RELAXED_ORBITALS", true);

        /*- change cavity mode polarization by redefining x, y, and z -*/
        options.add_str("ROTATE_POLARIZATION_AXIS", "XYZ");

        /*- residual norm -*/
        options.add_double("RESIDUAL_NORM",1.0e-5);

        /*- initial size of Davidson subspace (will be multiplied by number of desired roots) -*/
        options.add_int("INDIM", 5);

        /*- maximum size of Davidson subspace (will be multiplied by number of desired roots) -*/
        options.add_int("MAXDIM", 20);

        /*- number of roots -*/
        options.add_int("NUMBER_ROOTS", 5);


        /*- SUBSECTION CC_Cavity -*/

        /*- Maximum number of vectors in DIIS -*/
        options.add_int("DIIS_MAX_VECS", 8);

        /*- tile size for tiled array -*/
        options.add_int("TILE_SIZE", -1); // -1 means to use all data in one tile

        /*- number of threads to use for MADNESS -*/
        options.add_int("MAD_NUM_THREADS", 1);

        /*- do include triples in polaritioinic cc? -*/
        options.add_bool("QED_CC_INCLUDE_T3", false);
        options.add_bool("QED_CC_INCLUDE_U3", false);

        /*- do include quadruples in polaritioinic cc? -*/
        options.add_bool("QED_CC_INCLUDE_T4", false);
        options.add_bool("QED_CC_INCLUDE_U4", false);

        /*- do include u0 in polaritioinic ccsd? -*/
        options.add_bool("QED_CC_INCLUDE_U0", false);

        /*- do include u1 in polaritioinic ccsd? -*/
        options.add_bool("QED_CC_INCLUDE_U1", false);

        /*- do include u2 in polaritioinic ccsd? -*/
        options.add_bool("QED_CC_INCLUDE_U2", false);

        /*- explicitly set u0=0 -*/
        options.add_bool("ZERO_U0", false);

        /*- Perform EOM for cc_cavity -*/
        options.add_bool("PERFORM_EOM", false);

        /*- Perform EOM with full hamiltonian for qed-ccsd (not implemented) -*/
        options.add_bool("BUILD_HAMILTONIAN", false);

        /*- Convergence criteria for davidson solver -*/
        options.add_double("EOM_R_CONV", 1e-6);
        options.add_double("EOM_E_CONV", 1e-8);
        options.add_int("EOM_MAXITER", 100);

        /*- Use diagonals from singles block in EOM-CC -*/
        options.add_bool("EOM_SS_GUESS", true);

        /*- value of shift to apply to davidson solver -*/
        options.add_double("EOM_SHIFT", 0.0l); // 0.0 means no shift

        /*- Save and load eigenvectors -*/
        options.add_bool("SAVE_EVECS", false); // save eigenvectors to file?
        options.add_int("LOAD_ID", -1); // id of file to load eigenvectors from

        /*- EOM-CC TYPE -*/
        options.add_str("EOM_TYPE", "EE", "EE EA");

        /*- do print excited state transition dipoles (computed either way) -*/
        options.add_bool("EXCITED_PROPERTIES", true);

        /*- Select left and right states to visualize density or transition density -*/
        options.add_array("RDM_STATES");
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

    }else if ( options.get_str("HILBERT_METHOD") == "JELLIUM_SCF") {

        std::shared_ptr<Jellium_SCFSolver> jellium (new Jellium_SCFSolver(options));
        double energy = jellium->compute_energy();
        return ref_wfn;

    }else if ( options.get_str("HILBERT_METHOD") == "PP2RDM") {

        std::shared_ptr<pp2RDMSolver> pp2rdm (new pp2RDMSolver(ref_wfn,options));
        double energy = pp2rdm->compute_energy();
        return (std::shared_ptr<Wavefunction>)pp2rdm;

    }else if ( options.get_str("HILBERT_METHOD") == "P2RDM") {

        std::shared_ptr<p2RDMSolver> p2rdm (new p2RDMSolver(ref_wfn,options));
        double energy = p2rdm->compute_energy();
        return (std::shared_ptr<Wavefunction>)p2rdm;

    }else if ( options.get_str("HILBERT_METHOD") == "V2RDM_DOCI") {

        std::shared_ptr<v2RDM_DOCISolver> v2rdm_doci (new v2RDM_DOCISolver(ref_wfn,options));
        double energy = v2rdm_doci->compute_energy();
        return (std::shared_ptr<Wavefunction>)v2rdm_doci;

    }else if ( options.get_str("HILBERT_METHOD") == "V2RDM_CASSCF") {

        std::shared_ptr<v2RDMSolver> v2rdm (new v2RDMSolver(ref_wfn,options));
        double energy = v2rdm->compute_energy();

        if ( options.get_str("DERTYPE") == "FIRST" ) {

            // backtransform the tpdm
            std::vector<std::shared_ptr<MOSpace> > spaces;
            spaces.push_back(MOSpace::all);
            std::shared_ptr<TPDMBackTransform> transform = std::shared_ptr<TPDMBackTransform>(
            new TPDMBackTransform(ref_wfn,
                            spaces,
                            IntegralTransform::TransformationType::Unrestricted, // Transformation type
                            IntegralTransform::OutputType::DPDOnly,              // Output buffer
                            IntegralTransform::MOOrdering::QTOrder,              // MO ordering
                            IntegralTransform::FrozenOrbitals::None));           // Frozen orbitals?
            transform->backtransform_density();
            transform.reset();
        }

        return (std::shared_ptr<Wavefunction>)v2rdm;

    }else if ( options.get_str("HILBERT_METHOD") == "POLARITONIC_RHF") {

        std::shared_ptr<PolaritonicRHF> rhf (new PolaritonicRHF(ref_wfn,options));
        double energy = rhf->compute_energy();

        return (std::shared_ptr<Wavefunction>)rhf;

    }else if ( options.get_str("HILBERT_METHOD") == "POLARITONIC_TDDFT") {

        std::shared_ptr<PolaritonicRKS> rks (new PolaritonicRKS(ref_wfn,options));
        double energy = rks->compute_energy();

        std::shared_ptr<PolaritonicRTDDFT> rtddft (new PolaritonicRTDDFT((std::shared_ptr<Wavefunction>)rks,options,ref_wfn));
        double dum = rtddft->compute_energy();

        return (std::shared_ptr<Wavefunction>)rks;

    }else if ( options.get_str("HILBERT_METHOD") == "POLARITONIC_RCIS") {

        std::shared_ptr<PolaritonicRHF> rhf (new PolaritonicRHF(ref_wfn,options));
        double energy = rhf->compute_energy();

        std::shared_ptr<PolaritonicRCIS> rcis (new PolaritonicRCIS((std::shared_ptr<Wavefunction>)rhf,options));
        double dum = rcis->compute_energy();

        return (std::shared_ptr<Wavefunction>)rhf;

    }else if ( options.get_str("HILBERT_METHOD") == "POLARITONIC_UHF") {

        std::shared_ptr<PolaritonicUHF> uhf (new PolaritonicUHF(ref_wfn,options));
        double energy = uhf->compute_energy();

        return (std::shared_ptr<Wavefunction>)uhf;

    }else if ( options.get_str("HILBERT_METHOD") == "POLARITONIC_ROHF") {

        std::shared_ptr<PolaritonicROHF> rohf (new PolaritonicROHF(ref_wfn,options));
        double energy = rohf->compute_energy();

        return (std::shared_ptr<Wavefunction>)rohf;

    }else if ( options.get_str("HILBERT_METHOD") == "POLARITONIC_RKS") {

        std::shared_ptr<PolaritonicRKS> rks (new PolaritonicRKS(ref_wfn,options));
        double energy = rks->compute_energy();

        return (std::shared_ptr<Wavefunction>)rks;


    }else if ( options.get_str("HILBERT_METHOD") == "POLARITONIC_UKS") {

        std::shared_ptr<PolaritonicUKS> uks (new PolaritonicUKS(ref_wfn,options));
        double energy = uks->compute_energy();

        return (std::shared_ptr<Wavefunction>)uks;

    }else if ( options.get_str("HILBERT_METHOD") == "POLARITONIC_UCCSD") {

        if ( options.get_str("REFERENCE") == "UHF") {

            std::shared_ptr<PolaritonicUHF> uhf (new PolaritonicUHF(ref_wfn,options));
            double energy = uhf->compute_energy();

            std::shared_ptr<PolaritonicUCCSD> uccsd (new PolaritonicUCCSD((std::shared_ptr<Wavefunction>)uhf,options));
            energy = uccsd->compute_energy();
            return (std::shared_ptr<Wavefunction>)uccsd;

        }else if ( options.get_str("REFERENCE") == "ROHF" ) {

            std::shared_ptr<PolaritonicROHF> rohf (new PolaritonicROHF(ref_wfn,options));
            double energy = rohf->compute_energy();

            std::shared_ptr<PolaritonicUCCSD> uccsd (new PolaritonicUCCSD((std::shared_ptr<Wavefunction>)rohf,options));
            energy = uccsd->compute_energy();
            return (std::shared_ptr<Wavefunction>)uccsd;

        }else if ( options.get_str("REFERENCE") == "RHF" ) {

            std::shared_ptr<PolaritonicRHF> rhf (new PolaritonicRHF(ref_wfn,options));
            double energy = rhf->compute_energy();

            std::shared_ptr<PolaritonicUCCSD> uccsd (new PolaritonicUCCSD((std::shared_ptr<Wavefunction>)rhf,options));
            energy = uccsd->compute_energy();
            return (std::shared_ptr<Wavefunction>)uccsd;

        }else {

            throw PsiException("unknown REFERENCE for polaritonic UHF",__FILE__,__LINE__);

        }

    } else if ( options.get_str("HILBERT_METHOD") == "CC_CAVITY") {
#ifndef USE_QED_CC
        throw PsiException("CC_CAVITY requires the USE_QED_CC flag to be set at compile time",__FILE__,__LINE__);
#else

        // get the reference wavefunction
        std::shared_ptr<Wavefunction> qed_ref_wfn;
        if ( options.get_str("REFERENCE") == "UHF") {

            std::shared_ptr<PolaritonicUHF> uhf (new PolaritonicUHF(ref_wfn,options));
            double energy = uhf->compute_energy();

            qed_ref_wfn = (std::shared_ptr<Wavefunction>)uhf;

        } else if ( options.get_str("REFERENCE") == "ROHF" ) {

            std::shared_ptr<PolaritonicROHF> rohf (new PolaritonicROHF(ref_wfn,options));
            double energy = rohf->compute_energy();

            qed_ref_wfn = (std::shared_ptr<Wavefunction>)rohf;

        } else if ( options.get_str("REFERENCE") == "RHF" ) {

            std::shared_ptr<PolaritonicRHF> rhf (new PolaritonicRHF(ref_wfn,options));
            double energy = rhf->compute_energy();

            qed_ref_wfn = (std::shared_ptr<Wavefunction>)rhf;

        } else {
            throw PsiException("unknown REFERENCE for polaritonic UHF",__FILE__,__LINE__);
        }

        // initialize the tiled array library
        int argc = 1;
        char **argv = nullptr;
        psi::Process::environment.globals["MAD_NUM_THREADS"] = options.get_int("MAD_NUM_THREADS");

        // initialize tiledattay with MPIComm from mpi4py
        TA::initialize(argc, argv, CavityHelper::comm_);

        // create the CC_CAVITY object
        std::shared_ptr<CC_Cavity> qedcc;

        // determine the level of theory
        bool three_body = options.get_bool("QED_CC_INCLUDE_T3") || options.get_bool("QED_CC_INCLUDE_U3");
        bool four_body  = options.get_bool("QED_CC_INCLUDE_T4") || options.get_bool("QED_CC_INCLUDE_U4");
        bool has_photon = options.get_int("N_PHOTON_STATES") > 1;

        // select the appropriate derived CC_CAVITY object
        if (four_body)
    #if MAX_CC_LEVEL >= 4
            qedcc = std::shared_ptr<CC_Cavity>(new QED_CCSDTQ(qed_ref_wfn, options));
    #else
            throw PsiException("QED_CCSDTQ requires MAX_CC_LEVEL >= 4 at compile time",__FILE__,__LINE__);
    #endif
        else if (three_body)
    #if MAX_CC_LEVEL >= 3
            qedcc = std::shared_ptr<CC_Cavity>(new QED_CCSDT(qed_ref_wfn, options));
    #else
            throw PsiException("QED_CCSDT requires MAX_CC_LEVEL >= 3 at compile time",__FILE__,__LINE__);
    #endif
        else {
            if (has_photon)
                qedcc = std::shared_ptr<CC_Cavity>(new QED_CCSD(qed_ref_wfn, options));
            else
                qedcc = std::shared_ptr<CC_Cavity>(new CC_Cavity(qed_ref_wfn, options));
        }

        // compute the energy
        double energy = qedcc->compute_energy();

        // perform QED-EOM-CC if requested
        bool do_eom = options.get_bool("PERFORM_EOM");

        if (do_eom && !four_body && !three_body) {
        #ifndef USE_QED_EOM_CC
            throw PsiException("QED-EOM-CC calculations require the `USE_QED_EOM_CC` flag to be set at compile time",__FILE__,__LINE__);
        #else
            std::shared_ptr<EOM_Driver> eom_driver;
            if (options.get_str("EOM_TYPE") == "EE") { // use EOM for excitation energies
                if (has_photon)
                    eom_driver = std::shared_ptr<EOM_Driver>(
                        new EOM_EE_QED_CCSD((std::shared_ptr<CC_Cavity>) qedcc, options));
                else
                    eom_driver = std::shared_ptr<EOM_Driver>(
                        new EOM_EE_CCSD((std::shared_ptr<CC_Cavity>) qedcc, options));
            }
            else if (options.get_str("EOM_TYPE") == "EA") { // use EOM for electron attachment
                throw PsiException("QED-EOM-EA-CC is not fuctional for release yet.", __FILE__, __LINE__);
                eom_driver = std::shared_ptr<EOM_Driver>(
                        new EOM_EA_Driver((std::shared_ptr<CC_Cavity>) qedcc, options));
            }
            else throw PsiException("EOM_TYPE not recognized. Please choose EE or EA.", __FILE__, __LINE__);

            // compute the excited state energies for the given EOM type
            eom_driver->compute_eom_energy();

            // compute the RDMs and oscillator strengths for the given EOM type
            if (options.get_str("EOM_TYPE") == "EE") {
                std::shared_ptr<EOM_EE_RDM> rdm(new EOM_EE_RDM((std::shared_ptr<EOM_Driver>) eom_driver, options));
                rdm->compute_eom_rdms();
                rdm->compute_oscillators();
                rdm->print_oscillators();
                if (options.exists("RDM_STATES")) {
                    vector<int> rdm_states = options.get_int_vector("RDM_STATES");
                    if (rdm_states.size() != 2 && rdm_states.size() != 1) {
                        throw PsiException("RDM_STATES must be a vector of length 1 or 2: provided size is " + std::to_string(rdm_states.size()), __FILE__, __LINE__);
                    }

                    for (auto state : rdm_states) {
                        if (state < 0 || state > eom_driver->M_) {
                            throw PsiException("RDM_STATES must be between 0 and the number of states.", __FILE__, __LINE__);
                        }
                    }

                    rdm->molden(rdm_states);
                }

            }
        #endif

        }

        // Tear down the parallel runtime_ for last energy call when FINALIZE = true
        TA::finalize();

        // return the wavefunction
        return (std::shared_ptr<Wavefunction>)qedcc;
#endif
    } else {
        throw PsiException("unknown HILBERT_METHOD",__FILE__,__LINE__);
    }

    return ref_wfn;
}

} // End namespaces

