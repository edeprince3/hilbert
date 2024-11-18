message("")


option(USE_QED_CC "Build QED-CC functionality" OFF)
option(USE_QED_EOM_CC "Build QED-EOM-CC functionality" OFF)
option(USE_QED_2RDM "Build 2RDMs for QED-EOM-CC functionality" OFF)
option(KEEP_NO_QED "Build separate equations for CC calculations without QED" OFF)
option(MAX_CC_LEVEL "CC level for QED-CC calculations" 2)
option(MAX_PHOTON_LEVEL "Photon level for QED-CC calculations" 2)

# if USE_QED_EOM_CC is set, USE_QED_CC must also be set
if (USE_QED_EOM_CC)
    set(USE_QED_CC ON)
endif ()


# add files for performing QED-CC calculations
set(qed_cc "")
set(qed_eom_cc "")
set(qed_eom_builds "")
set(qed_resid_builds "")

if (USE_QED_CC)
    message(STATUS "QED-CC will be built")
    add_definitions(-DUSE_QED_CC=1)
else ()
    message(STATUS "QED-CC and QED-EOM-CC will NOT be built")
    return()
endif ()

if (NOT MAX_CC_LEVEL OR MAX_CC_LEVEL LESS 2 OR MAX_CC_LEVEL GREATER 4)
    message(FATAL_ERROR "MAX_CC_LEVEL must be set to 2, 3, or 4")
    exit()
endif ()

if (NOT MAX_PHOTON_LEVEL OR MAX_PHOTON_LEVEL GREATER 2 OR MAX_PHOTON_LEVEL LESS 1)
    message(FATAL_ERROR "MAX_PHOTON_LEVEL must be set to 1 or 2")
    exit()
endif ()

message(STATUS "MAX_CC_LEVEL = ${MAX_CC_LEVEL}")
message(STATUS "MAX_PHOTON_LEVEL = ${MAX_PHOTON_LEVEL}")
add_definitions(-DMAX_CC_LEVEL=${MAX_CC_LEVEL})
add_definitions(-DMAX_PHOTON_LEVEL=${MAX_PHOTON_LEVEL})


if (USE_QED_EOM_CC)
    message(STATUS "QED-EOM-CC will be built")
    add_definitions(-DUSE_QED_EOM_CC=1)
else ()
    message(STATUS "QED-EOM-CC will NOT be built")
endif ()

if (USE_QED_2RDM)
    message(STATUS "QED-EOM-RDMs-CC will be built")
    add_definitions(-BUILD_QED_2RDM=1)
else ()
    message(STATUS "2RDMs for QED-EOM-CC will NOT be built")
endif ()

if (KEEP_NO_QED)
    message(STATUS "Separate equations for CC calculations without QED will be built")
    add_definitions(-DKEEP_NO_QED=1)
else ()
    message(STATUS "Separate equations for CC calculations without QED will NOT be built")
endif ()

message("")

# disable var tracking since it is slow
add_compile_options(-fno-var-tracking)

# fetch cmake files for tiledarray from ValeevGroup/kit-cmake
include(FetchContent)

# enable python bindings for tiledarray (requires shared library)
#set(TA_PYTHON ON)
#set(BUILD_SHARED_LIBS ON)

# get mpi components
message(STATUS "Finding MPI")
find_package(MPI REQUIRED COMPONENTS C CXX)

find_package(mpi4py)

# check if mpi4py was found
if (NOT mpi4py_FOUND)
    message(FATAL_ERROR "mpi4py not found; please install mpi4py with the command 'pip install mpi4py' or 'conda install mpi4py'")
endif ()

# set include directory for mpi4py
include_directories(${mpi4py_INCLUDE_DIRS})

set(BUILD_TESTING OFF CACHE BOOL "Build tests" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared library" FORCE)

# add tiledarray library
FetchContent_Declare(tiledarray
        GIT_REPOSITORY https://github.com/ValeevGroup/tiledarray.git
        GIT_TAG 248f85ce93e35d7e5326f35daceb4c262987c167
)

FetchContent_MakeAvailable(tiledarray)

#set_target_properties(tiledarray PROPERTIES CMAKE_BUILD_TYPE Release)
#set_target_properties(tiledarray PROPERTIES CMAKE_POSITION_INDEPENDENT_CODE ON)


# files for QED-CC
set(qed_cc
        src/cc_cavity/src/cc_cavity.cc
        src/cc_cavity/src/ccsd/ccsd.cc
        src/cc_cavity/src/qed_ccsd_21/qed_ccsd_21.cc

        src/cc_cavity/misc/ta_diis.cc
        src/cc_cavity/misc/timer.cc
        src/cc_cavity/misc/ta_helper.cc
)

# files to build the QED-CC equations (residuals)
set(qed_resid_builds
        src/cc_cavity/src/qed_ccsd_21/residuals/resid_21_1.cc
        src/cc_cavity/src/qed_ccsd_21/residuals/resid_21_2.cc
        src/cc_cavity/src/qed_ccsd_21/residuals/resid_21_3.cc
        src/cc_cavity/src/qed_ccsd_21/residuals/resid_21_4.cc
)
if (KEEP_NO_QED)
    set(qed_resid_builds
            ${qed_resid_builds}
            src/cc_cavity/src/ccsd/residuals/resid_00_1.cc
            src/cc_cavity/src/ccsd/residuals/resid_00_2.cc
            src/cc_cavity/src/ccsd/residuals/resid_00_3.cc
            src/cc_cavity/src/ccsd/residuals/resid_00_4.cc
    )
endif ()

if (MAX_PHOTON_LEVEL GREATER 1)
    set(qed_cc
            ${qed_cc}
            src/cc_cavity/src/qed_ccsd_22/qed_ccsd_22.cc
    )
    set(qed_resid_builds
            ${qed_resid_builds}
            src/cc_cavity/src/qed_ccsd_22/residuals/resid_22_1.cc
            src/cc_cavity/src/qed_ccsd_22/residuals/resid_22_2.cc
            src/cc_cavity/src/qed_ccsd_22/residuals/resid_22_3.cc
            src/cc_cavity/src/qed_ccsd_22/residuals/resid_22_4.cc
            src/cc_cavity/src/qed_ccsd_22/residuals/resid_22_5.cc
            src/cc_cavity/src/qed_ccsd_22/residuals/resid_22_6.cc
    )
endif ()


if (MAX_CC_LEVEL GREATER 2)
    set(qed_cc
            ${qed_cc}
    )
endif ()

if (MAX_CC_LEVEL GREATER 3)
    set(qed_cc
            ${qed_cc}
    )
endif ()

if (MAX_CC_LEVEL GREATER 4)
    message(FATAL_ERROR "MAX_CC_LEVEL must be set to 2, 3, or 4")
endif ()

# the files in residuals take a long time to compile with -fvar-tracking-assignments 
# so we disable it for these files
SET_SOURCE_FILES_PROPERTIES(${qed_resid_builds} PROPERTIES COMPILE_FLAGS -fno-var-tracking)

# add files for equation-of-motion QED-CC calculations
if (USE_QED_EOM_CC)
    set(qed_eom_cc
            src/cc_cavity/src/eom_driver.cc
            src/cc_cavity/src/eom_rdm.cc

            src/cc_cavity/src/ccsd/eom_ee_ccsd.cc
            src/cc_cavity/src/ccsd/eom_ea_ccsd.cc
            src/cc_cavity/src/ccsd/eom_ee_rdm.cc

            src/cc_cavity/src/qed_ccsd_21/eom_ee_qed_ccsd_21.cc
            src/cc_cavity/src/qed_ccsd_21/eom_ee_qed_rdm_21.cc
            src/cc_cavity/src/qed_ccsd_21/eom_ea_qed_ccsd_21.cc

    )
    set(qed_eom_builds
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
    )

    if (USE_QED_2RDM)
        set(qed_eom_builds
                ${qed_eom_builds}
                src/cc_cavity/src/qed_ccsd_21/rdms/rdm2_21_1.cc
                src/cc_cavity/src/qed_ccsd_21/rdms/rdm2_21_2.cc
                src/cc_cavity/src/qed_ccsd_21/rdms/rdm2_21_3.cc
                src/cc_cavity/src/qed_ccsd_21/rdms/rdm2_21_4.cc
                src/cc_cavity/src/qed_ccsd_21/rdms/rdm2_21_5.cc
                src/cc_cavity/src/qed_ccsd_21/rdms/rdm2_21_6.cc
                src/cc_cavity/src/qed_ccsd_21/rdms/rdm2_21_7.cc
                src/cc_cavity/src/qed_ccsd_21/rdms/rdm2_21_8.cc
        )
    endif ()

    if (KEEP_NO_QED)
        set(qed_eom_builds
                ${qed_eom_builds}
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

        )

        if (USE_QED_2RDM)
            set(qed_eom_builds
                    ${qed_eom_builds}
                    src/cc_cavity/src/ccsd/rdms/rdm2_00_1.cc
                    src/cc_cavity/src/ccsd/rdms/rdm2_00_2.cc
                    src/cc_cavity/src/ccsd/rdms/rdm2_00_3.cc
                    src/cc_cavity/src/ccsd/rdms/rdm2_00_4.cc
            )
        endif ()
    endif ()

    # the files in sigma_builds, eom_hamiltonians, and rdms take a long time to compile with -fvar-tracking-assignments
    # so we disable it for these files
    SET_SOURCE_FILES_PROPERTIES(${qed_eom_builds} PROPERTIES COMPILE_FLAGS -fno-var-tracking)


endif ()

# add the files to the qed_cc list
set(qed_cc ${qed_cc} ${qed_resid_builds} ${qed_eom_cc} ${qed_eom_builds})

message("")
