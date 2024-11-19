option(WITH_TA "Build QED-CC and EOM-QED-CC functionality with TiledArray" OFF)

# add files for performing QED-CC calculations
set(qed_cc "")
set(qed_eom_cc "")
set(qed_eom_builds "")
set(qed_resid_builds "")

if (WITH_TA)
    message(STATUS "QED-CC and EOM-QED-CC with TiledArray will be built")
    add_definitions(-DWITH_TA=ON)
else ()
    message(STATUS "QED-CC and EOM-QED-CC with TiledArray will NOT be built")
    return()
endif ()

# fetch cmake files for tiledarray from ValeevGroup/kit-cmake
include(FetchContent)

# get mpi components
message(STATUS "Finding MPI")
find_package(MPI REQUIRED COMPONENTS C CXX)

find_package(mpi4py)

# check if mpi4py was found
if (NOT mpi4py_FOUND)
    message(FATAL_ERROR "mpi4py not found; please install mpi4py with the command 'conda install mpi4py' or 'pip install mpi4py'")
endif ()

# add tiledarray library
message(STATUS "Fetching TiledArray")
FetchContent_Declare(tiledarray
        GIT_REPOSITORY https://github.com/ValeevGroup/tiledarray.git
        GIT_TAG 248f85ce93e35d7e5326f35daceb4c262987c167
)

# disable building the tests for tiledarray
set(BUILD_TESTING OFF CACHE BOOL "Build tests" FORCE)
FetchContent_MakeAvailable(tiledarray)

# Tiled Array should always be built in Release mode
set_target_properties(tiledarray PROPERTIES CMAKE_BUILD_TYPE Release)

# files for QED-CC
set(qed_cc
        src/cc_cavity/src/cc_cavity.cc
        src/cc_cavity/src/ccsd/ccsd.cc
        src/cc_cavity/src/qed_ccsd_21/qed_ccsd_21.cc
        src/cc_cavity/src/qed_ccsd_22/qed_ccsd_22.cc

        src/cc_cavity/misc/ta_diis.cc
        src/cc_cavity/misc/timer.cc
        src/cc_cavity/misc/ta_helper.cc
)

# files to build the QED-CC equations (residuals)
set(qed_resid_builds
        ${qed_resid_builds} # CCSD
        src/cc_cavity/src/ccsd/residuals/resid_00_1.cc
        src/cc_cavity/src/ccsd/residuals/resid_00_2.cc
        src/cc_cavity/src/ccsd/residuals/resid_00_3.cc
        src/cc_cavity/src/ccsd/residuals/resid_00_4.cc
)
set(qed_resid_builds # QED-CCSD-21
        ${qed_resid_builds}
        src/cc_cavity/src/qed_ccsd_21/residuals/resid_21_1.cc
        src/cc_cavity/src/qed_ccsd_21/residuals/resid_21_2.cc
        src/cc_cavity/src/qed_ccsd_21/residuals/resid_21_3.cc
        src/cc_cavity/src/qed_ccsd_21/residuals/resid_21_4.cc
)
set(qed_resid_builds
        ${qed_resid_builds} # QED-CCSD-22
        src/cc_cavity/src/qed_ccsd_22/residuals/resid_22_1.cc
        src/cc_cavity/src/qed_ccsd_22/residuals/resid_22_2.cc
        src/cc_cavity/src/qed_ccsd_22/residuals/resid_22_3.cc
        src/cc_cavity/src/qed_ccsd_22/residuals/resid_22_4.cc
        src/cc_cavity/src/qed_ccsd_22/residuals/resid_22_5.cc
        src/cc_cavity/src/qed_ccsd_22/residuals/resid_22_6.cc
)

# add files for equation-of-motion QED-CC calculations
set(qed_eom_cc
        src/cc_cavity/src/eom_driver.cc
        src/cc_cavity/src/eom_rdm.cc

        src/cc_cavity/src/ccsd/eom_ee_ccsd.cc
        src/cc_cavity/src/ccsd/eom_ea_ccsd.cc

        src/cc_cavity/src/qed_ccsd_21/eom_ee_qed_ccsd_21.cc
        src/cc_cavity/src/qed_ccsd_21/eom_ea_qed_ccsd_21.cc

)

set(qed_eom_builds
        ${qed_eom_builds} # CCSD
        src/cc_cavity/src/ccsd/sigma_builds/eom_ee/sigma_00_1.cc
        src/cc_cavity/src/ccsd/sigma_builds/eom_ee/sigma_00_2.cc
        src/cc_cavity/src/ccsd/sigma_builds/eom_ee/sigma_00_3.cc
        src/cc_cavity/src/ccsd/sigma_builds/eom_ee/sigma_00_4.cc
        src/cc_cavity/src/ccsd/sigma_builds/eom_ee/intermediates_00.cc

        src/cc_cavity/src/ccsd/sigma_builds/eom_ea/sigma_00_1.cc
        src/cc_cavity/src/ccsd/sigma_builds/eom_ea/sigma_00_2.cc
        src/cc_cavity/src/ccsd/sigma_builds/eom_ea/sigma_00_3.cc
        src/cc_cavity/src/ccsd/sigma_builds/eom_ea/sigma_00_4.cc
        src/cc_cavity/src/ccsd/sigma_builds/eom_ea/intermediates_00.cc

        src/cc_cavity/src/ccsd/eom_ee_rdm.cc
)
set(qed_eom_builds
        ${qed_eom_builds} # QED-CCSD-21
        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ee/intermediates_21.cc
        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ee/sigma_21_1.cc
        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ee/sigma_21_2.cc
        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ee/sigma_21_3.cc
        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ee/sigma_21_4.cc
        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ee/sigma_21_5.cc
        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ee/sigma_21_6.cc
        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ee/sigma_21_7.cc
        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ee/sigma_21_8.cc

        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ea/sigma_21_1.cc
        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ea/sigma_21_2.cc
        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ea/sigma_21_3.cc
        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ea/sigma_21_4.cc
        src/cc_cavity/src/qed_ccsd_21/sigma_builds/eom_ea/intermediates_21.cc

        src/cc_cavity/src/qed_ccsd_21/eom_ee_qed_rdm_21.cc
)

# the following files take a long time to compile with -fvar-tracking-assignments; disable it for these files
SET_SOURCE_FILES_PROPERTIES(${qed_resid_builds} PROPERTIES COMPILE_FLAGS -fno-var-tracking)
SET_SOURCE_FILES_PROPERTIES(${qed_eom_builds} PROPERTIES COMPILE_FLAGS -fno-var-tracking)

# add the files to the qed_cc list
set(qed_cc ${qed_cc} ${qed_resid_builds} ${qed_eom_cc} ${qed_eom_builds})
message("")