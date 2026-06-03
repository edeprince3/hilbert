# Building Hilbert on macOS Apple Silicon

This guide describes a native arm64 build of Hilbert as a Psi4 plugin on
Apple Silicon. It focuses on the QED-CC/TiledArray path, so it enables
`WITH_TA=ON` and uses MPI through the active conda environment.

Keep the whole stack native arm64: Python, Psi4, compilers, MPI,
BLAS/LAPACK, and the final `hilbert.so`. Avoid mixing conda MPI, Homebrew
MPI, and Rosetta/x86_64 tools in the same build.

## 1. Activate The Conda Environment

Use a native Apple Silicon terminal.

```bash
source /path/to/miniforge3/etc/profile.d/conda.sh
conda activate p4dev
```

Confirm that Python and the shell are arm64.

```bash
uname -m
python -c "import platform; print(platform.machine())"
```

Both commands should print `arm64`.

## 2. Install Build Dependencies

Install the build tools and libraries into the same environment that provides
Psi4.

```bash
conda install -c conda-forge cmake ninja compilers openmpi mpi4py pybind11 eigen boost-cpp pylibxc
```

The QED-CC/TiledArray build uses MPI and `mpi4py`. The MC-PDFT tests use
`pylibxc`.

## 3. Check Psi4 And MPI

Make sure `psi4`, Python, and MPI all come from the intended environment.

```bash
which python
which psi4
python --version
psi4 --version
python -c "import psi4; print(psi4.__file__)"
psi4 --plugin-compile

which mpicc
which mpicxx
which mpifort
python -c "import mpi4py; print('mpi4py OK')"
```

If you also have a local Psi4 source build, check `PYTHONPATH`. A path such as
`/path/to/psi4/objdir-Release/stage/lib` can shadow the conda Psi4 package. For
clean conda-based checks, unset it:

```bash
env -u PYTHONPATH python -c "import psi4; print(psi4.__file__)"
```

## 4. Configure

From the Hilbert repository root, use Psi4's plugin CMake cache. The cache path
is usually under `$CONDA_PREFIX/share/cmake/psi4`.

```bash
cmake \
  -C "$CONDA_PREFIX/share/cmake/psi4/psi4PluginCache.cmake" \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DCMAKE_INSTALL_PREFIX="$PWD" \
  -DWITH_TA=ON \
  -DFETCHCONTENT_UPDATES_DISCONNECTED=ON \
  -DCMAKE_C_COMPILER="$(which mpicc)" \
  -DCMAKE_CXX_COMPILER="$(which mpicxx)" \
  -DCMAKE_Fortran_COMPILER="$(which mpifort)" \
  -DCMAKE_CXX_FLAGS="-march=native -stdlib=libc++ -D_LIBCPP_ENABLE_CXX20_REMOVED_SHARED_PTR_UNIQUE" \
  -G Ninja \
  -S . \
  -B objdir
```

Notes:

- `-stdlib=libc++` keeps the plugin on the macOS C++ runtime used by the Psi4
  plugin build.
- `-D_LIBCPP_ENABLE_CXX20_REMOVED_SHARED_PTR_UNIQUE` is needed because MADNESS
  still uses `std::shared_ptr::unique()` while this Psi4 plugin build uses a
  C++20-style libc++ configuration.
- `-march=native` is fine for building and running on the same Mac. Omit it if
  the binary needs to run on a different Apple Silicon model.
- A fresh configure may need network access to populate TiledArray, MADNESS,
  Umpire, and related dependencies. If dependency sources have not been fetched
  yet and configure fails, rerun with network access and either omit
  `FETCHCONTENT_UPDATES_DISCONNECTED` or set it to `OFF`.

On macOS, CMake often selects Apple Accelerate for BLAS/LAPACK. That is fine
for a first build as long as the rest of the stack is consistently arm64.

## 5. Build And Install

```bash
cmake --build objdir -j
cmake --install objdir
```

For a readable first failure, rerun with one build job:

```bash
cmake --build objdir -j1
```

The install step places `hilbert.so` at the repository root when
`CMAKE_INSTALL_PREFIX="$PWD"` is used. It also installs dependency headers and
CMake metadata under root-level directories such as `include/`, `lib/`, and
`_deps/`; these are local build artifacts and are ignored by `.gitignore`.

## 6. Verify

Check the built plugin architecture.

```bash
file hilbert.so
```

Expected output includes `Mach-O 64-bit bundle arm64`.

Run a basic import check from outside the repository, making sure conda Psi4 is
used.

```bash
cd /tmp
env -u PYTHONPATH python - <<'PY'
import sys
import psi4
print(psi4.__file__)
sys.path.insert(0, "/path/to/parent/of/hilbert")
import hilbert
print(hilbert.__version__)
PY
```

Then run a quick functional test.

```bash
cd /path/to/hilbert/tests
env -u PYTHONPATH make quick
```

## Source Compatibility Notes

Two small source edits were needed for this macOS/Psi4 1.10-style build:

- `src/cc_cavity/misc/timer.cc` now includes `<sstream>` explicitly before
  using `std::stringstream`. Newer libc++ headers do not guarantee this through
  transitive includes.
- `src/cc_cavity/misc/ta_helper.cc` now calls `array.init_elements(...)`
  directly. The previous `array.template init_elements(...)` form is invalid
  without an explicit template argument list and is rejected by current Clang in
  this C++20-style plugin build.
