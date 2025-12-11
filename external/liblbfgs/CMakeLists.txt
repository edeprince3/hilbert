
find_package(LBFGS CONFIG QUIET)

include(FetchContent)

message(STATUS "liblbfgs will be built at compile time")
FetchContent_Declare(liblbfgs_external
    GIT_REPOSITORY https://github.com/edeprince3/liblbfgs.git

)
FetchContent_MakeAvailable(liblbfgs_external)

set(liblbfgs_DIR ${EXTERNAL_INSTALL_PREFIX}/share/cmake/liblbfgs CACHE PATH "path to internally built liblbfgsConfig.cmake" FORCE)
