"""
Helper function for CQED_RHF

"""

__authors__   = ["Jon McTague", "Jonathan Foley"]
__credits__   = ["Jon McTague", "Jonathan Foley"]

__copyright_amp__ = "(c) 2014-2018, The Psi4NumPy Developers"
__license__   = "BSD-3-Clause"
__date__      = "2021-01-15"

# ==> Import Psi4, NumPy, & SciPy <==
import psi4
import numpy as np
import scipy.linalg as la
import time

def cqed_rhf(lambda_vector, molecule_string, self_consistent_dipole=True):
    """ Computes the QED-RHF energy and density 

        Arguments
        ---------
        lambda_vector : 1 x 3 array of floats
            the electric field vector
        
        molecule_string : string
            specifies the molecular geometry

        self_consistent_dipole : bool (optional)
            specifies if the dipole moment expectation value should be updated within the SCF iterations
            'False' will simply use the RHF dipole expectation value throughout
            defaults to True if no argument passed

        Returns
        -------
        psi4_rhf_energy : float
            Ground state energy from canonical RHF wavefunction from

        cqed_rhf_energy : float
            Ground state energy of the CQED_RHF wavefunction

        cqed_c : nmo x nmo array of floats
            Transformation vectors corresponding to the CQED_RHF orbitals

        Example
        -------
        >>> psi4_rhf_energy, cqed_rhf_energy, cqed_rhf_orbs = cqed_rhf([0., 0., 1e-2], '''\nMg\nH 1 1.7\nsymmetry c1\n1 1\n''')
        
    """
    # define geometry using the molecule_string
    mol = psi4.geometry(molecule_string)
    # run psi4 to get ordinary scf energy and wavefunction object
    psi4_rhf_energy, wfn = psi4.energy('scf', return_wfn=True)

    # Create instance of MintsHelper class
    mints = psi4.core.MintsHelper(wfn.basisset())

    # Grab data from wavfunction
    # number of doubly occupied orbitals
    ndocc   = wfn.nalpha()

    # grab all transformation vectors and store to a numpy array!
    C = np.asarray(wfn.Ca())

    # Compute required quantities for SCF
    V = np.asarray(mints.ao_potential())
    T = np.asarray(mints.ao_kinetic())
    I = np.asarray(mints.ao_eri())

    # Extra terms for Pauli-Fierz Hamiltonian
    # nuclear dipole
    mu_nuc_x = mol.nuclear_dipole()[0]
    mu_nuc_y = mol.nuclear_dipole()[1]
    mu_nuc_z = mol.nuclear_dipole()[2]

    # dipole arrays in AO basis
    mu_ao_x = np.asarray(mints.ao_dipole()[0])
    mu_ao_y = np.asarray(mints.ao_dipole()[1])
    mu_ao_z = np.asarray(mints.ao_dipole()[2])

    # transform dipole array to canonical MO basis from ordinary RHF (no photon)
    mu_cmo_x = np.dot(C.T, mu_ao_x).dot(C)
    mu_cmo_y = np.dot(C.T, mu_ao_y).dot(C)
    mu_cmo_z = np.dot(C.T, mu_ao_z).dot(C)

    # compute components of electronic dipole moment <mu> from ordinary RHF (no photon)
    mu_exp_x = 0.0
    mu_exp_y = 0.0
    mu_exp_z = 0.0
    for i in range(0, ndocc):
        # double because this is only alpha terms!
        mu_exp_x += 2 * mu_cmo_x[i, i]
        mu_exp_y += 2 * mu_cmo_y[i, i]
        mu_exp_z += 2 * mu_cmo_z[i, i]

    # need to add the nuclear term to the expectation values above which
    # only included the electronic term!
    mu_exp_x += mu_nuc_x
    mu_exp_y += mu_nuc_y
    mu_exp_z += mu_nuc_z

    # We need to carry around the electric field dotted into the nuclear dipole moment
    # and the electric field dotted into the RHF electronic dipole expectation value...
    # \lambda_vector \cdot \mu_{nuc}
    l_dot_mu_nuc = lambda_vector[0] * mu_nuc_x + lambda_vector[1] * mu_nuc_y + lambda_vector[2] * mu_nuc_z
    # \lambda_vecto \cdot < \mu > where <\mu> contains electronic and nuclear contributions
    l_dot_mu_exp = lambda_vector[0] * mu_exp_x + lambda_vector[1] * mu_exp_y + lambda_vector[2] * mu_exp_z


    # dipole constants to add to E_RHF
    #  0.5 * (\lambda_vector \cdot \mu_{nuc})** 2 
    #      - (\lambda_vector \cdot <\mu> ) ( \lambda_vector\cdot \mu_{nuc})
    # +0.5 * (\lambda_vector \cdot <\mu>) ** 2
    d_c = 0.5 * l_dot_mu_nuc **2 - l_dot_mu_nuc * l_dot_mu_exp + 0.5 * l_dot_mu_exp ** 2

    # quadrupole arrays
    Q_ao_xx = np.asarray(mints.ao_quadrupole()[0])
    Q_ao_xy = np.asarray(mints.ao_quadrupole()[1])
    Q_ao_xz = np.asarray(mints.ao_quadrupole()[2])
    Q_ao_yy = np.asarray(mints.ao_quadrupole()[3])
    Q_ao_yz = np.asarray(mints.ao_quadrupole()[4])
    Q_ao_zz = np.asarray(mints.ao_quadrupole()[5])

    # ordinary H_core
    H_0 = T + V

    # Pauli-Fierz 1-e quadrupole terms ... these terms have a factor of 1/2
    Q_PF = 0.5 * lambda_vector[0] * lambda_vector[0] * Q_ao_xx
    Q_PF += 0.5 * lambda_vector[1] * lambda_vector[1] * Q_ao_yy
    Q_PF += 0.5 * lambda_vector[2] * lambda_vector[2] * Q_ao_zz

    # accounting for the fact that Q_ij = Q_ji
    # by weighting Q_ij x 2... which cancels factor of 1/2
    Q_PF += lambda_vector[0] * lambda_vector[1] * Q_ao_xy
    Q_PF += lambda_vector[0] * lambda_vector[2] * Q_ao_xz
    Q_PF += lambda_vector[1] * lambda_vector[2] * Q_ao_yz

    # Pauli-Fierz 1-e dipole terms scaled by (\lambda_vector \cdot \mu_{nuc} - \lambda_vector \cdot <\mu>)
    d_PF =  (l_dot_mu_nuc - l_dot_mu_exp) * lambda_vector[0] * mu_ao_x
    d_PF += (l_dot_mu_nuc - l_dot_mu_exp) * lambda_vector[1] * mu_ao_y
    d_PF += (l_dot_mu_nuc - l_dot_mu_exp) * lambda_vector[2] * mu_ao_z

    # Add Pauli-Fierz terms to H_core
    H = H_0 + Q_PF + d_PF

    # Overlap for DIIS
    S = mints.ao_overlap()
    # Orthogonalizer A = S^(-1/2) using Psi4's matrix power.
    A = mints.ao_overlap()
    A.power(-0.5, 1.0e-16)
    A = np.asarray(A)

    # Calculate initial core guess: [Szabo:1996] pp. 145
    Hp = A.dot(H).dot(A)  # Eqn. 3.177
    e, C2 = np.linalg.eigh(Hp)  # Solving Eqn. 1.178
    C = A.dot(C2)  # Back transform, Eqn. 3.174
    Cocc = C[:, :ndocc]

    D = np.einsum("pi,qi->pq", Cocc, Cocc)  # [Szabo:1996] Eqn. 3.145, pp. 139

    #print("\nTotal time taken for setup: %.3f seconds" % (time.time() - t))

    print("\nStart SCF iterations:\n")
    t = time.time()
    E = 0.0
    Enuc = mol.nuclear_repulsion_energy()
    Eold = 0.0
    Dold = np.zeros_like(D)

    E_1el = np.einsum("pq,pq->", H + H, D) + Enuc + d_c
    print("One-electron energy = %4.16f" % E_1el)

    # Set convergence criteria
    maxiter = 200
    E_conv = 1.0e-6
    D_conv = 5.0e-5
    t = time.time()
    for SCF_ITER in range(1, maxiter + 1):

        # Build fock matrix: [Szabo:1996] Eqn. 3.154, pp. 141
        J = np.einsum("pqrs,rs->pq", I, D)
        K = np.einsum("prqs,rs->pq", I, D)

        # Pauli-Fierz dipole-dipole matrices
        M_xx = np.einsum("pq,rs,rs->pq", lambda_vector[0] * mu_ao_x, lambda_vector[0] * mu_ao_x, D)
        M_yy = np.einsum("pq,rs,rs->pq", lambda_vector[1] * mu_ao_y, lambda_vector[1] * mu_ao_y, D)
        M_zz = np.einsum("pq,rs,rs->pq", lambda_vector[2] * mu_ao_z, lambda_vector[2] * mu_ao_z, D)

        M_xy = np.einsum("pq,rs,rs->pq", lambda_vector[0] * mu_ao_x, lambda_vector[1] * mu_ao_y, D)
        M_xz = np.einsum("pq,rs,rs->pq", lambda_vector[0] * mu_ao_x, lambda_vector[2] * mu_ao_z, D)
        M_yz = np.einsum("pq,rs,rs->pq", lambda_vector[1] * mu_ao_y, lambda_vector[2] * mu_ao_z, D)

        # Pauli-Fierz dipole-dipole "exchange" terms
        N_xx = np.einsum("pr,qs,rs->pq", lambda_vector[0] * mu_ao_x, lambda_vector[0] * mu_ao_x, D)
        N_yy = np.einsum("pr,qs,rs->pq", lambda_vector[1] * mu_ao_y, lambda_vector[1] * mu_ao_y, D)
        N_zz = np.einsum("pr,qs,rs->pq", lambda_vector[2] * mu_ao_z, lambda_vector[2] * mu_ao_z, D)

        N_xy = np.einsum("pr,qs,rs->pq", lambda_vector[0] * mu_ao_x, lambda_vector[1] * mu_ao_y, D)
        N_xz = np.einsum("pr,qs,rs->pq", lambda_vector[0] * mu_ao_x, lambda_vector[2] * mu_ao_z, D)
        N_yz = np.einsum("pr,qs,rs->pq", lambda_vector[1] * mu_ao_y, lambda_vector[2] * mu_ao_z, D)

        # Build fock matrix: [Szabo:1996] Eqn. 3.154, pp. 141 
        # 
        # Pauli-Fierz terms
        F = H + J * 2 - K
        F += M_xx
        F += M_yy
        F +=  M_zz

        F += 2 * M_xy
        F += 2 * M_xz
        F += 2 * M_yz

        F -= 0.5 * N_xx
        F -= 0.5 * N_yy
        F -= 0.5 * N_zz

        F -= N_xy
        F -= N_xz
        F -= N_yz
        
        
        diis_e = np.einsum("ij,jk,kl->il", F, D, S) - np.einsum("ij,jk,kl->il", S, D, F)
        diis_e = A.dot(diis_e).dot(A)
        dRMS = np.mean(diis_e ** 2) ** 0.5

        # SCF energy and update: [Szabo:1996], Eqn. 3.184, pp. 150
        SCF_E = np.einsum("pq,pq->", F + H, D) + Enuc + d_c

        #print(
        #    "SCF Iteration %3d: Energy = %4.16f   dE = % 1.5E   dRMS = %1.5E"
        #    % (SCF_ITER, SCF_E, (SCF_E - Eold), dRMS)
        #)
        if (abs(SCF_E - Eold) < E_conv) and (dRMS < D_conv):
            break

        Eold = SCF_E
        Dold = D

        # Diagonalize Fock matrix: [Szabo:1996] pp. 145
        Fp = A.dot(F).dot(A)  # Eqn. 3.177
        e, C2 = np.linalg.eigh(Fp)  # Solving Eqn. 1.178
        C = A.dot(C2)  # Back transform, Eqn. 3.174
        Cocc = C[:, :ndocc]
        D = np.einsum("pi,qi->pq", Cocc, Cocc)  # [Szabo:1996] Eqn. 3.145, pp. 139

        if self_consistent_dipole:
            # transform dipole array to current MO basis
            mu_cmo_x = np.dot(C.T, mu_ao_x).dot(C)
            mu_cmo_y = np.dot(C.T, mu_ao_y).dot(C)
            mu_cmo_z = np.dot(C.T, mu_ao_z).dot(C)
               
            # compute components of electronic dipole moment <mu> 
            mu_exp_x = 0.0
            mu_exp_y = 0.0
            mu_exp_z = 0.0
               
            for i in range(0, ndocc):
                # double because this is only alpha terms!
                mu_exp_x += 2 * mu_cmo_x[i, i]
                mu_exp_y += 2 * mu_cmo_y[i, i]
                mu_exp_z += 2 * mu_cmo_z[i, i]
            # update <\mu> in current CQED-RHF basis 
            mu_exp_x += mu_nuc_x
            mu_exp_y += mu_nuc_y
            mu_exp_z += mu_nuc_z
                
            # update \lambda \cdot <\mu>
            l_dot_mu_exp = lambda_vector[0] * mu_exp_x + lambda_vector[1] * mu_exp_y + lambda_vector[2] * mu_exp_z
            d_PF =  (l_dot_mu_nuc - l_dot_mu_exp) * lambda_vector[0] * mu_ao_x
            d_PF += (l_dot_mu_nuc - l_dot_mu_exp) * lambda_vector[1] * mu_ao_y
            d_PF += (l_dot_mu_nuc - l_dot_mu_exp) * lambda_vector[2] * mu_ao_z

            # update Core Hamiltonian
            H = H_0 + Q_PF + d_PF

            # update dipole energetic contribution
            d_c = 0.5 * l_dot_mu_nuc **2 - l_dot_mu_nuc * l_dot_mu_exp + 0.5 * l_dot_mu_exp ** 2
        

        if SCF_ITER == maxiter:
            psi4.core.clean()
            raise Exception("Maximum number of SCF cycles exceeded.")
            
    print("Performed QED-RHF on the following molecule")
    print(molecule_string)
    print("Total time for SCF iterations: %.3f seconds \n" % (time.time() - t))

    print("QED-RHF   energy: %.8f hartree" % SCF_E)
    print("Psi4  SCF energy: %.8f hartree" % psi4_rhf_energy)

    return psi4_rhf_energy, SCF_E, C
