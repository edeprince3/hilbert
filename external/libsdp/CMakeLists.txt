
include(FetchContent)

message(STATUS "libsdp will be built at compile time")
FetchContent_Declare(libsdp_external
    GIT_REPOSITORY https://github.com/edeprince3/libsdp.git
    GIT_TAG 56a755f7ba63b8f14e06db81ee2d2409bf81c95c
)
FetchContent_MakeAvailable(libsdp_external)

set(libsdp_DIR ${EXTERNAL_INSTALL_PREFIX}/share/cmake/libsdp CACHE PATH "path to internally built libsdpConfig.cmake" FORCE)

