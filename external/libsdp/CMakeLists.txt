
include(FetchContent)

message(STATUS "libsdp will be built at compile time")
FetchContent_Declare(libsdp_external
    GIT_REPOSITORY https://github.com/edeprince3/libsdp.git
    GIT_TAG 132ecccb197abb3199a71f5fd6cead1f0102a6df
)
FetchContent_MakeAvailable(libsdp_external)

set(libsdp_DIR ${EXTERNAL_INSTALL_PREFIX}/share/cmake/libsdp CACHE PATH "path to internally built libsdpConfig.cmake" FORCE)

