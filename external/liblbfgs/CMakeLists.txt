
find_package(LBFGS CONFIG QUIET)

if(${LBFGS_FOUND})
    get_property(_loc TARGET LBFGS::lib PROPERTY LOCATION)
    message(STATUS "${Cyan}Found LBFGS${ColourReset}: ${_loc} (found version ${LBFGS_VERSION})")
    add_library(lbfgs_external INTERFACE)  # dummy
else()
    include(ExternalProject)
    message(STATUS "liblbfgs will be built at compile time")
    ExternalProject_Add(liblbfgs_external
        GIT_REPOSITORY https://github.com/chokkan/liblbfgs.git
    
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}

        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    )
    
    set(liblbfgs_DIR ${EXTERNAL_INSTALL_PREFIX}/share/cmake/liblbfgs CACHE PATH "path to internally built liblbfgsConfig.cmake" FORCE)
endif()

