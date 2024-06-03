message("") 


option(USE_QED_CC "Build QED-CC functionality" OFF)
option(USE_QED_EOM_CC "Build QED-EOM-CC functionality" OFF)

# if USE_QED_EOM_CC is set, USE_QED_CC must also be set
if (USE_QED_EOM_CC)
    set(USE_QED_CC ON)
endif()


# add files for performing QED-CC calculations
set(qed_cc "")
if (NOT USE_QED_CC)
  message(STATUS "QED-CC and QED-EOM-CC will NOT be built")
  return()
else()
  message(STATUS "QED-CC will be built")
  add_definitions(-DUSE_QED_CC=1)
endif()

if (USE_QED_EOM_CC)
  message(STATUS "QED-EOM-CC will be built")
  add_definitions(-DUSE_QED_EOM_CC=1)
else()
  message(STATUS "QED-EOM-CC will NOT be built")
endif()

# disable var tracking since it is slow
add_compile_options(-fno-var-tracking)

# add option for MAX_CC_LEVEL (default is 2, but can be set to 3 or 4)
option(MAX_CC_LEVEL "CC level for QED-CC calculations" 2)

# fetch cmake files for tiledarray from ValeevGroup/kit-cmake
include(FetchContent)
FetchContent_Declare(kit-cmake
  GIT_REPOSITORY https://github.com/ValeevGroup/kit-cmake.git
)
FetchContent_MakeAvailable(kit-cmake)

# include cmake modules and toolchains from kit-cmake
list(APPEND CMAKE_MODULE_PATH 
    ${kit-cmake_SOURCE_DIR}/modules
    ${kit-cmake_SOURCE_DIR}/toolchains
)

# use toolchain file from kit-cmake depending on the compiler
option(CUSTOM_TOOLCHAIN "name of custom toolchain file from kit-cmake" "")

if (CUSTOM_TOOLCHAIN)
  include(${kit-cmake_SOURCE_DIR}/toolchains/${CUSTOM_TOOLCHAIN})
else()
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
      message(STATUS "Using GCC toolchain")  
      include(${kit-cmake_SOURCE_DIR}/toolchains/gcc-mpi-mkl-tbb.cmake)
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        message(STATUS "Using Intel toolchain")
        include(${kit-cmake_SOURCE_DIR}/toolchains/intel-parallel-studio.cmake)
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
      message(STATUS "Using Intel oneAPI toolchain")
      include(${kit-cmake_SOURCE_DIR}/toolchains/intel-oneapi)
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(STATUS "Using Clang toolchain")
        include(${kit-cmake_SOURCE_DIR}/toolchains/clang-mpi-mkl-tbb.cmake)
    endif()
endif()

# set initial compile flags
set(CMAKE_CXX_FLAGS_INIT           "" CACHE STRING "Initial C++ compile flags")
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g -Wall" CACHE STRING "Initial C++ debug compile flags")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -mcpu=native -DNDEBUG" CACHE STRING "Initial C++ minimum size release compile flags")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -mcpu=native -DNDEBUG" CACHE STRING "Initial C++ release compile flags")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -Wall" CACHE STRING "Initial C++ release with debug info compile flags")

# print build type
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
  
# enable python bindings for tiledarray (requires shared library)
#set(TA_PYTHON ON)
#set(BUILD_SHARED_LIBS ON)

# get mpi components
find_package(MPI REQUIRED)

# get include directory for mpi4py
find_package(MPI4PY REQUIRED)

# add tiledarray library
FetchContent_Declare(tiledarray
  GIT_REPOSITORY https://github.com/ValeevGroup/tiledarray.git
  GIT_TAG b2254f352260d80e15c5f4a94b7bfcad71f118fb
)
FetchContent_MakeAvailable(tiledarray)

set_target_properties(tiledarray PROPERTIES CMAKE_BUILD_TYPE Release)
set_target_properties(tiledarray PROPERTIES CMAKE_POSITION_INDEPENDENT_CODE ON)

# files for QED-CC
set(qed_cc
  src/cc_cavity/src/cc_cavity.cc
  src/cc_cavity/src/derived/qed_ccsd.cc
  src/cc_cavity/src/derived/lambda_driver.cc

  src/cc_cavity/misc/ta_diis.cc
  src/cc_cavity/misc/timer.cc
  src/cc_cavity/misc/ta_helper.cc
)

# files to build the QED-CC equations (residuals)
set(qed_cc_builds
  src/cc_cavity/src/derived/residuals/ccsd_resid.cc
  src/cc_cavity/src/derived/residuals/qed_ccsd_resid.cc
  src/cc_cavity/src/derived/residuals/qed_ccsd_lam_resid.cc
)

# the files in residuals take a long time to compile with -fvar-tracking-assignments 
# so we disable it for these files
SET_SOURCE_FILES_PROPERTIES( ${qed_cc_builds} PROPERTIES COMPILE_FLAGS -fno-var-tracking)

set(qed_cc ${qed_cc} ${qed_cc_builds})


if (NOT MAX_CC_LEVEL OR MAX_CC_LEVEL LESS 2)
  set(MAX_CC_LEVEL 2)
endif()

if (MAX_CC_LEVEL GREATER 2)
  set(qed_cc
    ${qed_cc}
    # src/cc_cavity/src/derived/qed_ccsdt.cc
    # src/cc_cavity/src/derived/residuals/qed_ccsdt_resid.cc
  )
endif()

if (MAX_CC_LEVEL GREATER 3)
  set(qed_cc
    ${qed_cc}
    # src/cc_cavity/src/derived/qed_ccsdtq.cc
    # src/cc_cavity/src/derived/residuals/qed_ccsdtq_resid.cc
  )
endif()

add_definitions(-DMAX_CC_LEVEL=${MAX_CC_LEVEL})

# add files for equation-of-motion QED-CC calculations
if (USE_QED_EOM_CC)
  set(qed_eom_cc
    src/cc_cavity/src/eom_driver.cc
    src/cc_cavity/src/eom_rdm.cc
    src/cc_cavity/src/derived/eom_ee_ccsd.cc
    src/cc_cavity/src/derived/eom_ee_qed_ccsd.cc
    src/cc_cavity/src/derived/eom_ea_driver.cc
    src/cc_cavity/src/derived/eom_ee_rdm.cc
  )
  set(qed_eom_builds
    src/cc_cavity/src/derived/sigma_builds/eom_ee_qed_ccsd_build.cc
    src/cc_cavity/src/derived/sigma_builds/eom_ee_qed_ccsd_intermediates.cc
    src/cc_cavity/src/derived/sigma_builds/eom_ee_ccsd_build.cc
    src/cc_cavity/src/derived/sigma_builds/eom_ea_build.cc

    src/cc_cavity/src/derived/eom_hamiltonians/eom_ea_hamiltonian.cc

    src/cc_cavity/src/derived/2rdms/eom_ee_2rdm.cc
  )

  # the files in sigma_builds, eom_hamiltonians, and rdms take a long time to compile with -fvar-tracking-assignments 
  # so we disable it for these files
  SET_SOURCE_FILES_PROPERTIES( ${qed_eom_builds} PROPERTIES COMPILE_FLAGS -fno-var-tracking)

  # add the files to the qed_cc list
  set(qed_cc ${qed_cc} ${qed_eom_cc} ${qed_eom_builds})

endif()
message("")
