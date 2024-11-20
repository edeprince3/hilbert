# [CC-CAVITY](#CC-CAVITY)


## Introduction
The `cc_cavity` class is a module of Hilbert which implements the Coupled-Cluster (CC) and 
Equation-of-Motion Coupled-Cluster (EOM-CC) methods with quantum electrodynamics (QED) for the simulation of 
cavity polaritons using the [TiledArray](https://github.com/ValeevGroup/tiledarrayhttps://github.com/ValeevGroup/tiledarray) library.

## Features

- Implements QED-CC methods for ground state calculations
  - CCSD and QED-CCSD-1 methods are implemented
  - CCSDT, CCSDTQ, QED-CCSDT-1, and QED-CCSDTQ-1 infrastructures are in place, however, the residual equations are not yet open sourced.


- Provides QED-EOM-CC methods for excited state calculations (optional)
  - EOM-EE-CCSD and QED-EOM-EE-CCSD-1 methods are implemented (excitation energies)
  - EOM-EA-CCSD and QED-EOM-EA-CCSD-1 infrastructures are in place, however, it is still under development.
  - The EOM-CC methods are implemented using the Non-Symmetric Davidson algorithm, which is an interative subspace diagonalization method.
- Builds the 1-RDMs for the excited state methods, which is used to compute the transition dipole moments.
  - Only available for EOM-EE-CCSD and QED-EOM-EE-CCSD-1 methods.
  - 2-RDMs are not yet implemented.

- Utilizes the TiledArray library for efficient tensor operations
  - Supports MPI distributed memory parallelization, however, this is not yet implemented for EOM-CC methods.
- Includes ex

## Directory Structure

- `include`: Contains header files for the `cc_cavity` class
  - `derived`: Contains header files for derived classes
- `src`: Contains the source files for the `cc_cavity` class, including derived classes for QED-CC and QED-EOM-CC methods
    - `derived`: Contains source files for derived classes
- `misc`: Contains miscellaneous helper functions and utilities

## Installation

Follow the same instructions as the README.md file in the root directory of the Hilbert project, 
with the following modifications:

- To enable QED-CC functionality, set the `WITH_TA` flag in CMake to `ON`.
  ``` bash
  cmake {...psi4-generated stuff...} -DWITH_TA=ON -Bobjdir -DCMAKE_INSTALL_PREFIX=path_to_hilbert_top_dir
  ```
- it is also recommended to add `-G Ninja` to the cmake line to speed up the build process. 
  It is also more stable for installing TiledArray, which is a dependency of CC_Cavity.

After completing these steps, the `cc_cavity` class will be available as a plugin in your Psi4 installation. 
You can then use it to perform QED-CC and QED-EOM-CC calculations.

## Dependencies

- TiledArray, and its [dependencies](https://github.com/ValeevGroup/tiledarray/blob/master/INSTALL.md#prerequisites) (automatically built by CMake if not available)
- Psi4 (which this plugin is built for) and its [dependencies](https://psicode.org/psi4manual/master/build_planning.html#faq-coredepend)
- MPI (required even with a single node)
- MPI4PY
  - This is a Python package that is used to interface with the MPI library through Psi4.
  - This package must be installed using the same python interpreter that Psi4 is using.
  - If you are using the Psi4 conda environment, then you can install this package using the following command:
  `conda install -c conda-forge mpi4py`
  - If not using conda, you can check which python interpreter Psi4 uses by running the following command:
    ``` bash
    grep -m 1 pyexe_dir `which psi4` | awk -F'"' '{print $2}'
    ```
  - Install the mpi4py package using the same python interpreter:
  `/path/to/python -m pip install mpi4py`

### Usage

After installing the `cc_cavity` plugin, you can use it within Psi4 to perform QED-CC and QED-EOM-CC calculations. To do so, follow these steps:

1. Prepare your input file for Psi4, including the molecular geometry, basis set, and other desired options. 
2. Add `import hilbert` to the top of your input file and add the path to your Hilbert installation directory to the `PYTHONPATH` environment variable or using the `sys` module in your input file.
3. Set the desired method and options for the `cc_cavity` calculation using the `set` command.
4. Run the calculation with `energy('cc_cavity')`.

### Example
Here is an example input file to run with Psi4 that has comments for the different settings that can be included.

``` python
sys.path.insert(0, '../')
import hilbert

molecule {
  0 1
  O
  H 1 0.96
  H 1 0.96 2 104.5
  
  # do not reorient the molecule.
  # The orientation relative to the cavity polarization affects the result
  no_reorient 
  no_com
  
  # symmetry is not supported for QED-CC calculations
  symmetry c1
}

set {
  basis cc-pvdz # basis set for the molecular orbitals
  reference uhf # all calculations use a spin-orbital basis
  
  # CD algorithm is recommended for QED-CC calculations. DF is also supported, but not PK. 
  scf_type cd 
  cholesky_tolerance 1e-12
  
  # convergence criteria for the SCF and CC iterations
  e_convergence 1e-10
  r_convergence 1e-10
  d_convergence 1e-10       
}
 
lam = 0.05 # coupling factor (this value is for a 0.74 nm³ cavity)
frequency = 1.0 # frequency of the cavity mode
coupling_strength = (lam / math.sqrt(2.0 * frequency)) # coupling strength
set hilbert {
  maxiter 1000 # maximum number of iterations for the groundstate CC methods
  
  # set the cavity parameters
  cavity_frequency         [0.0, 0.0, $frequency]
  cavity_coupling_strength [0.0, 0.0, $coupling_strength]
  
  QED_USE_RELAXED_ORBITALS  true # do relax the SCF orbitals within the cavity
  ROTATE_POLARIZATION_AXIS   XYZ # rotate the cavity polarization from the XYZ axis to ZXY, or YZX.
  
  NUMBER_ROOTS                20 # number of excited states to compute
  MAXDIM                       8 # maximum dimension of the EOM-CC subspace (multiplied by NUMBER_ROOTS)
  INDIM                        4 # initial dimension of the EOM-CC subspace (multiplied by NUMBER_ROOTS)
  
  EOM_MAXITER                250 # maximum number of iterations for the EOM-CC iterations
  MAD_NUM_THREADS              1 # the number of MADNESS threads to use with mpirun 
                                 # Note: if not using mpi or using eom, set this to 1 for best performanc
  TILE_SIZE                   -1 # the spacing between tiles in a tiledarray (-1 places all data on a single tile)
                                 # Note: if not using mpi or using eom, set this to -1 for best performance
}

# memory for Psi4 (note: this is for the integrals, not the CC calculation which does not restrict memory usage)
memory 1000 MB

# run the QED-CC calculations

en1, wfn = energy('qed-ccsd', return_wfn=True)        # QED-CCSD-21 w/o TiledArray
en2, wfn = energy('qed-ccsd-21', return_wfn=True)     # QED-CCSD-21 w/ TiledArray
en3, wfn = energy('qed-ccsd-22', return_wfn=True)     # QED-CCSD-22 w/ TiledArray

# print the energies
print('QED-CCSD-21: ', en1)
print('QED-CCSD-21 (TiledArray): ', en2)
print('QED-CCSD-22 (TiledArray): ', en3)

set hilbert EOM_TYPE EE
en4, wfn = energy('eom-qed-ccsd-21', return_wfn=True) # EOM-EE-QED-CCSD-21 w/ TiledArray

set hilbert EOM_TYPE EA
en5, wfn = energy('eom-qed-ccsd-21', return_wfn=True) # EOM-EA-QED-CCSD-21 w/ TiledArray
```

If the name of this file is `input.dat`, the file can be run with the following command:
``` bash
psi4 -n $NUM_THREADS input.dat output.dat
```

Where $NUMTHREADS is the number of threads you wish to run the calculation with (when using parallel BLAS).
To run the calculation with mpi, you can use the following command:
``` bash
mpirun -n $NUM_THREADS psi4 input.dat output.dat
```
Note: MPI functionality is not optimized as of now.