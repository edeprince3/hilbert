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

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(psi4OptionsTools)

find_package(pybind11 REQUIRED)
find_package(psi4 1.3 REQUIRED)

add_psi4_plugin(hilbert

    src/plugin.cc

    src/python_api/python_helpers.cc

    src/doci/david.cc
    src/doci/doci_solver.cc
    src/doci/write_tpdm.cc

    src/pp2rdm/minimization_routines.cc
    src/pp2rdm/pp2rdm_solver.cc
    src/pp2rdm/write_tpdm.cc

    src/misc/threeindexintegrals.cc
    src/misc/cg_solver.cc
    src/misc/diis.cc

    src/focas/update_transformation_matrix.cc
    src/focas/orbital_optimizer.cc

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

    src/v2rdm_doci/backtransform_tpdm.cc
    src/v2rdm_doci/basis.cc
    src/v2rdm_doci/checkpoint.cc
    src/v2rdm_doci/d2.cc
    src/v2rdm_doci/d3.cc
    src/v2rdm_doci/g2.cc
    src/v2rdm_doci/integraltransform_sort_so_tpdm.cc
    src/v2rdm_doci/integraltransform_tpdm_unrestricted.cc
    src/v2rdm_doci/oei.cc
    src/v2rdm_doci/orbital_lagrangian.cc
    src/v2rdm_doci/q2.cc
    src/v2rdm_doci/sortintegrals.cc
    src/v2rdm_doci/t1.cc
    src/v2rdm_doci/t2.cc
    src/v2rdm_doci/tei.cc
    src/v2rdm_doci/update_primal.cc
    src/v2rdm_doci/v2rdm_doci_solver.cc
    src/v2rdm_doci/write_tpdm.cc
    src/v2rdm_doci/write_tpdm_iwl.cc

    src/v2rdm_casscf/backtransform_tpdm.cc
    src/v2rdm_casscf/basis.cc
    src/v2rdm_casscf/checkpoint.cc
    src/v2rdm_casscf/d2.cc
    src/v2rdm_casscf/d3.cc
    src/v2rdm_casscf/d4.cc
    src/v2rdm_casscf/extended_koopmans.cc
    src/v2rdm_casscf/fcidump.cc
    src/v2rdm_casscf/g2.cc
    src/v2rdm_casscf/integraltransform_sort_so_tpdm.cc
    src/v2rdm_casscf/integraltransform_tpdm_unrestricted.cc
    src/v2rdm_casscf/natural_orbitals.cc
    src/v2rdm_casscf/oei.cc
    src/v2rdm_casscf/orbital_lagrangian.cc
    #src/v2rdm_casscf/python_api.cc
    src/v2rdm_casscf/q2.cc
    src/v2rdm_casscf/sortintegrals.cc
    src/v2rdm_casscf/spin.cc
    src/v2rdm_casscf/t1.cc
    src/v2rdm_casscf/t2.cc
    src/v2rdm_casscf/tei.cc
    src/v2rdm_casscf/update_primal.cc
    #src/v2rdm_casscf/v2rdm_helper.cc
    src/v2rdm_casscf/v2rdm_solver.cc
    src/v2rdm_casscf/write_3pdm.cc
    src/v2rdm_casscf/write_opdm.cc
    src/v2rdm_casscf/write_tpdm.cc
    src/v2rdm_casscf/write_tpdm_iwl.cc
)

# external projects 
set(EXTERNAL_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/external/)
add_subdirectory(external)

# link liblbfgs and find header
set(LBFGS_LIBRARY "${PROJECT_SOURCE_DIR}/external/lib/liblbfgs.a")
target_link_libraries(hilbert PRIVATE ${LBFGS_LIBRARY})
set(EXTERNAL_INCLUDE "${PROJECT_SOURCE_DIR}/external/include")
target_include_directories(hilbert PRIVATE ${EXTERNAL_INCLUDE})

# hilbert module headers
target_include_directories(hilbert PRIVATE "${PROJECT_SOURCE_DIR}/src")

add_dependencies(hilbert liblbfgs_external)

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
