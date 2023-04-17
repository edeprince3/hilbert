if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()

cmake_policy(SET CMP0048 NEW)  # project_VERSION* variables populated from project(... VERSION x.x.x) string

project(hilbert
        VERSION 0.1
        LANGUAGES CXX C Fortran)

set(hilbert_AUTHORS      "A. Eugene DePrince III and the DePrince Group")
set(hilbert_DESCRIPTION  "Hilbert: A space for quantum chemistry plugins to Psi4")
set(hilbert_URL          "https://github.com/edeprince3/hilbert")
set(hilbert_LICENSE      "GPL-3.0+")
set(CMAKE_CXX_STANDARD 17)

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(pybind11 REQUIRED)
find_package(psi4 REQUIRED)
include(psi4OptionsTools)

# add files for performing EOM-CC calculations if requested
option(USE_QED_EOM_CC "Build QED-EOM-CC functionality" OFF)
option(USE_QED_CC "Build QED-CC functionality" OFF)

# add option for MAX_CC_LEVEL (default is 2, but can be set to 3 or 4)
set(MAX_CC_LEVEL 2 CACHE STRING "CC level for QED-CC calculations")

# add files for performing QED-CC calculations
set(qed_cc "")
set(qed_eom_cc "")

if (USE_QED_CC)
  message(STATUS "QED-CC will be built")
  if (USE_QED_EOM_CC)
    message(STATUS "QED-EOM-CC will be built")
  else()
    message(STATUS "QED-EOM-CC will NOT be built")
  endif()
  message(STATUS "MAX_CC_LEVEL is set to ${MAX_CC_LEVEL}")
else()
    message(STATUS "QED-CC will NOT be built")
    if (USE_QED_EOM_CC)
        message(FATAL_ERROR "USE_QED_EOM_CC requires USE_QED_CC")
    endif()
endif()

if (USE_QED_CC)

  # find prerequisites
  if(Fortran_COMPILER_ID MATCHES "Intel" AND NOT I_MPI_ROOT)
    # find_package(MPI) doesn't find mpi fortran compilers from intel and hangs if it's not set
    message(FATAL_ERROR "I_MPI_ROOT must be set to the root directory of the Intel MPI installation")
  endif()
  find_package(Boost REQUIRED COMPONENTS serialization)
  find_package(MPI REQUIRED)
  find_package(Eigen3 REQUIRED)
  find_package(OpenMP)

  # get include directory for mpi4py
  if (NOT MPI4PY_INCLUDE_DIR)
    if (NOT Python_EXECUTABLE) # We need to find Python if it's not already found
      find_package(Python COMPONENTS Interpreter REQUIRED)
    endif()
    execute_process( # Get the mpi4py include directory (this is a pip package and does not have a CMake config file)
            COMMAND
            ${Python_EXECUTABLE} -c "import mpi4py; import os; print(mpi4py.get_include(), end='')"
            OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
    )
    if(MPI4PY_INCLUDE_DIR STREQUAL "")
      message(FATAL_ERROR "Could not find mpi4py include directory. Install mpi4py with pip or set MPI4PY_INCLUDE_DIR.")
    endif()
  endif()

  # add include directories
  message(STATUS "Found mpi4py include directory: ${MPI4PY_INCLUDE_DIR}")
  include_directories(${MPI4PY_INCLUDE_DIR})

  # add tiledarray library
  include(FetchContent)
  FetchContent_Declare( tiledarray
          GIT_REPOSITORY https://github.com/ValeevGroup/tiledarray.git
          GIT_TAG c42361ce7f455eb3d77a86396d502f37c0b771e8
          )
  FetchContent_GetProperties(tiledarray)
  if(NOT tiledarray_POPULATED)
    FetchContent_Populate(tiledarray)
    add_subdirectory(${tiledarray_SOURCE_DIR} ${tiledarray_BINARY_DIR})
  endif()

  # add files for performing QED-CC calculations
  set(qed_cc
          src/cc_cavity/src/cc_cavity.cc
          src/cc_cavity/src/derived/qed_ccsd.cc
          src/cc_cavity/src/derived/qed_ccsd_resid.cc

          src/cc_cavity/misc/diis_qed.cc
          src/cc_cavity/misc/nonsym_davidson_solver_qed.cc
          src/cc_cavity/misc/threeindexintegralsta.cc
          src/cc_cavity/misc/timer.cc
          src/cc_cavity/misc/ta_helper.cc
          )
  add_definitions(-DUSE_QED_CC=1)

  # set default MAX_CC_LEVEL to 2 if it is not set
  if (NOT MAX_CC_LEVEL OR MAX_CC_LEVEL LESS 2)
    set(MAX_CC_LEVEL 2)
  endif()

  # add files for performing QED-CC calculations at MAX_CC_LEVEL >= 3
  if (MAX_CC_LEVEL GREATER 2)
      set(qed_cc
              ${qed_cc}
              src/cc_cavity/src/derived/qed_ccsdt.cc
              src/cc_cavity/src/derived/qed_ccsdt_resid.cc
              )
  endif()

  # add files for performing QED-CC calculations at MAX_CC_LEVEL >= 4
  if (MAX_CC_LEVEL GREATER 3)
      set(qed_cc
              ${qed_cc}
              src/cc_cavity/src/derived/qed_ccsdtq.cc
              src/cc_cavity/src/derived/qed_ccsdtq_resid.cc
              )
  endif()

  # add MAX_CC_LEVEL to compiler flags
  add_definitions(-DMAX_CC_LEVEL=${MAX_CC_LEVEL})
endif()

if (USE_QED_EOM_CC AND USE_QED_CC)
  set(qed_eom_cc
          src/cc_cavity/src/eom_driver.cc
          src/cc_cavity/src/eom_rdm.cc

          src/cc_cavity/src/derived/eom_ee_driver.cc
          src/cc_cavity/src/derived/eom_ee_build.cc
          src/cc_cavity/src/derived/eom_ea_driver.cc
          src/cc_cavity/src/derived/eom_ea_build.cc
          src/cc_cavity/src/derived/eom_ee_rdm.cc
          )
  add_definitions(-DUSE_QED_EOM_CC=1)
endif()

# files for performing polaritonic calculations
set(polaritonic_scf
        src/polaritonic_scf/hf.cc
        src/polaritonic_scf/rhf.cc
        src/polaritonic_scf/uhf.cc
        src/polaritonic_scf/rohf.cc
        src/polaritonic_scf/rks.cc
        src/polaritonic_scf/uks.cc
        src/polaritonic_scf/rcis.cc
        src/polaritonic_scf/rtddft.cc
        src/polaritonic_scf/utddft.cc
        src/polaritonic_scf/uccsd.cc
        )
set(jellium
        src/jellium/jellium_scf_solver.cc
        src/jellium/jellium_integrals.cc
        src/jellium/legendre.cc
        )
set(doci
        src/doci/doci_solver.cc
        src/doci/write_tpdm.cc
        src/doci/fcidump.cc
        )
set(p2rdm
        src/p2rdm/p2rdm_solver.cc
        src/p2rdm/p2rdm_residual.cc
        src/p2rdm/write_tpdm.cc
        )
set(pp2rdm
        src/pp2rdm/pp2rdm_solver.cc
        src/pp2rdm/minimization_routines.cc
        src/pp2rdm/write_tpdm.cc
        )
set(misc
        src/misc/real_space_density.cc
        src/misc/backtransform_tpdm.cc
        src/misc/davidson_solver.cc
        src/misc/nonsym_davidson_solver.cc
        src/misc/diis.cc
        src/misc/lbfgs_helper.cc
        src/misc/integraltransform_sort_so_tpdm.cc
        src/misc/integraltransform_tpdm_unrestricted.cc
        src/misc/threeindexintegrals.cc
        )
set(focas
        src/focas/update_transformation_matrix.cc
        src/focas/orbital_optimizer.cc
        src/focas/C_files/orbopt_energy.cc
        src/focas/C_files/orbopt_gradient.cc
        src/focas/C_files/orbopt_hessian.cc
        src/focas/C_files/orbopt_exponential.cc
        src/focas/C_files/orbopt_integraltransform_cpu.cc
        src/focas/C_files/orbopt_sort.cc
        src/focas/C_files/orbopt_optimize.cc
        src/focas/C_files/orbopt_utility.cc
        src/focas/C_files/orbopt_lagrangian.cc

        src/focas/focas_data.F90
        src/focas/focas_exponential.F90
        src/focas/focas_hessian.F90
        src/focas/focas_transform_driver.F90
        src/focas/focas_diis.F90
        src/focas/focas_genfock.F90
        src/focas/focas_interface.F90
        src/focas/focas_transform_oeints.F90
        src/focas/focas_driver.F90
        src/focas/focas_gradient.F90
        src/focas/focas_redundant.F90
        src/focas/focas_transform_teints.F90
        src/focas/focas_energy.F90
        src/focas/focas_gradient_hessian.F90
        src/focas/focas_semicanonical.F90
        src/focas/focas_full_hessian.F90
        )

set(v2rdm_doci
        src/v2rdm_doci/basis.cc
        src/v2rdm_doci/checkpoint.cc
        src/v2rdm_doci/d2.cc
        src/v2rdm_doci/d3.cc
        src/v2rdm_doci/g2.cc
        src/v2rdm_doci/oei.cc
        src/v2rdm_doci/orbital_lagrangian.cc
        src/v2rdm_doci/q2.cc
        src/v2rdm_doci/sortintegrals.cc
        src/v2rdm_doci/t1.cc
        src/v2rdm_doci/t2.cc
        src/v2rdm_doci/tei.cc
        src/v2rdm_doci/v2rdm_doci_solver.cc
        src/v2rdm_doci/write_tpdm.cc
        src/v2rdm_doci/write_tpdm_iwl.cc
        )

set(v2rdm_casscf
        src/v2rdm_casscf/basis.cc
        src/v2rdm_casscf/checkpoint.cc
        src/v2rdm_casscf/d2.cc
        src/v2rdm_casscf/d3.cc
        src/v2rdm_casscf/e3.cc
        src/v2rdm_casscf/f3.cc
        src/v2rdm_casscf/q3.cc
        src/v2rdm_casscf/d4.cc
        src/v2rdm_casscf/extended_koopmans.cc
        src/v2rdm_casscf/fcidump.cc
        src/v2rdm_casscf/g2.cc
        src/v2rdm_casscf/generalized_pauli.cc
        src/v2rdm_casscf/generalized_pauli_3_10.cc
        src/v2rdm_casscf/generalized_pauli_3_6.cc
        src/v2rdm_casscf/generalized_pauli_3_8.cc
        src/v2rdm_casscf/generalized_pauli_4_10.cc
        src/v2rdm_casscf/generalized_pauli_4_8.cc
        src/v2rdm_casscf/generalized_pauli_5_10.cc
        src/v2rdm_casscf/natural_orbitals.cc
        src/v2rdm_casscf/oei.cc
        src/v2rdm_casscf/orbital_lagrangian.cc
        src/v2rdm_casscf/python_api.cc
        src/v2rdm_casscf/q2.cc
        src/v2rdm_casscf/sortintegrals.cc
        src/v2rdm_casscf/spin.cc
        src/v2rdm_casscf/t1.cc
        src/v2rdm_casscf/t2.cc
        src/v2rdm_casscf/tei.cc
        src/v2rdm_casscf/update_transformation_matrix.cc
        src/v2rdm_casscf/v2rdm_solver.cc
        src/v2rdm_casscf/write_3pdm.cc
        src/v2rdm_casscf/write_opdm.cc
        src/v2rdm_casscf/write_tpdm.cc
        src/v2rdm_casscf/write_tpdm_iwl.cc
        )


add_psi4_plugin(hilbert

        src/plugin.cc
        src/python_api/python_helpers.cc

        ${qed_cc}
        ${qed_eom_cc}
        ${polaritonic_scf}
        ${jellium}
        ${doci}
        ${p2rdm}
        ${pp2rdm}
        ${misc}
        ${focas}
        ${v2rdm_doci}
        ${v2rdm_casscf}
        )

# external projects 
set(EXTERNAL_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/external/)
add_subdirectory(external)

set(EXTERNAL_INCLUDE  "${PROJECT_SOURCE_DIR}/external/include")
target_include_directories(hilbert PRIVATE  ${EXTERNAL_INCLUDE})
   
# hilbert module headers
target_include_directories(hilbert PRIVATE "${PROJECT_SOURCE_DIR}/src")

target_link_libraries(hilbert PUBLIC sdp)
target_link_libraries(hilbert PUBLIC LBFGS::lib)

if(USE_QED_CC)
  target_include_directories(hilbert PRIVATE ${MPI4PY_INCLUDE_DIR})
  target_include_directories(hilbert PRIVATE ${tiledarray_INCLUDE_DIR})
  target_link_libraries(hilbert PUBLIC MPI::MPI_CXX MPI::MPI_C MPI::MPI_Fortran tiledarray)
endif()

install(DIRECTORY ${hilbert_BINARY_DIR}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        USE_SOURCE_PERMISSIONS
        PATTERN "*.pyc" EXCLUDE
        PATTERN "*.mod" EXCLUDE
        PATTERN "CMakeCache.txt" EXCLUDE
        PATTERN "CMakeFiles" EXCLUDE
        PATTERN "cmake_install.cmake" EXCLUDE
        PATTERN "Makefile" EXCLUDE
        PATTERN "install_manifest.txt" EXCLUDE
)
