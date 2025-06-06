cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(hilbert
        VERSION 0.1
        LANGUAGES CXX C Fortran)

set(hilbert_AUTHORS      "A. Eugene DePrince III and the DePrince Group")
set(hilbert_DESCRIPTION  "Hilbert: A space for quantum chemistry plugins to Psi4")
set(hilbert_URL          "https://github.com/edeprince3/hilbert")
set(hilbert_LICENSE      "GPL-3.0+")
set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(psi4OptionsTools)

find_package(psi4 REQUIRED)
find_package(pybind11 REQUIRED)

# Include QED-CC CMake configurations if requested
include(${PROJECT_SOURCE_DIR}/cmake/cc_cavity.cmake)


#### groups for different calculation types ###

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
        src/pp2rdm/minimization_routines.cc
        src/pp2rdm/pp2rdm_solver.cc
        src/pp2rdm/write_tpdm.cc
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

set(misc
        src/misc/real_space_density.cc
        src/misc/backtransform_tpdm.cc
        src/misc/diis.cc
        src/misc/lbfgs_helper.cc
        src/misc/integraltransform_sort_so_tpdm.cc
        src/misc/integraltransform_tpdm_unrestricted.cc
        src/misc/threeindexintegrals.cc
        src/misc/schmidt_add.cc
        src/misc/sq_rsp.cc
        src/misc/davidson_solver.cc
        src/misc/nonsym_davidson_solver.cc
)

add_psi4_plugin(hilbert

        src/plugin.cc
        src/python_api/python_helpers.cc

        ${qed_cc}
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

if(WITH_TA)
  target_link_libraries(hilbert PUBLIC mpi4py::mpi4py tiledarray)
endif()

# set install prefix to project root in hilbert folder
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/hilbert)
    message(STATUS "Setting install prefix to ${CMAKE_INSTALL_PREFIX}")
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
