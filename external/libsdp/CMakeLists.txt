
#find_package(SDP CONFIG QUIET)

#if(${SDP})
#    get_property(_loc TARGET SDP::lib PROPERTY LOCATION)
#    message(STATUS "${Cyan}Found SDP${ColourReset}: ${_loc} (found version ${SDP_VERSION})")
#    add_library(sdp_external INTERFACE)  # dummy
#else()
    include(ExternalProject)
    message(STATUS "libsdp will be built at compile time")
    #set(CATCH_GIT_TAG 5a2e83dc4a58bd191a5c4f5145e507bf0a03ba38)
    # this one
    #set(CATCH_GIT_TAG 00b5ff68f7e7d0914eeab6b1de668cb0f1e6fc81)
    #set(CATCH_GIT_TAG 933ba6227f534b9d7940ccd16420940ce1ee24d0)
    set(CATCH_GIT_TAG 24047db4503565ce9b10140abfd7d29f9ad78bc8) # most recent as of 10/26/21
    ExternalProject_Add(libsdp_external
        GIT_REPOSITORY https://github.com/edeprince3/libsdp.git
        GIT_TAG ${CATCH_GIT_TAG}
    
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                   -DLBFGS_INCLUDE=${EXTERNAL_INSTALL_PREFIX}/include

        CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    )
    add_dependencies(libsdp_external liblbfgs_external)

    set(libsdp_DIR ${EXTERNAL_INSTALL_PREFIX}/share/cmake/libsdp CACHE PATH "path to internally built libsdpConfig.cmake" FORCE)
#endif()

