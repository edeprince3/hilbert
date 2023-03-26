
include(FetchContent)

message(STATUS "libsdp will be built at compile time")
FetchContent_Declare(libsdp_external
    GIT_REPOSITORY https://github.com/edeprince3/libsdp.git
    #GIT_TAG 9b54d6fbb894cb48971dcbb1f4e41b55ea1f8f4c
    GIT_TAG 33803b906a327f0c27bb329c32812820fdc6c57c
)
FetchContent_MakeAvailable(libsdp_external)

set(libsdp_DIR ${EXTERNAL_INSTALL_PREFIX}/share/cmake/libsdp CACHE PATH "path to internally built libsdpConfig.cmake" FORCE)

