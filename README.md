# OSLO — Unified ODE Solver Toolbox

![License](https://img.shields.io/github/license/MarcoGrossi92/OSLO)
![Fortran](https://img.shields.io/badge/Fortran-90%2B-blue)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey)

**OSLO (ODE SoLver toolbOx)** is a high-performance, Fortran-based toolbox for solving systems of ordinary differential equations (ODEs).
It provides a **unified and modular interface** to multiple state-of-the-art ODE solvers, allowing users to switch between stiff and non-stiff integrators with minimal code changes.

OSLO is designed for **scientific computing and HPC environments**, bundling well-established solvers libraries under a single, consistent API.

---

## Key Features

* **Unified Solver Interface**
  Switch between different ODE solvers using a common Fortran API.

* **Broad Solver Coverage**

  * Hairer & Wanner solvers for stiff and non-stiff problems
  * **FATODE** with automatic differentiation
  * **SUNDIALS** (CVODE, IDA, etc.)

* **Parallel Execution**

  * Shared-memory parallelism via **OpenMP**
  * Distributed-memory parallelism via **MPI**

* **Compiler Flexibility**

  * GNU Fortran (`gfortran`)
  * Intel Fortran (`ifort` / `ifx`)

* **Robust Linear Algebra**

  * LAPACK integration for dense and sparse systems

* **Modular Architecture**

  * Easy to extend or embed in existing Fortran codes

---

## Requirements

### Mandatory

* **CMake** ≥ 3.23
* **Fortran compiler**: GNU or Intel
* **LAPACK**
* **C compiler** (required for SUNDIALS)

### Optional

* **OpenMP** — shared-memory parallelism
* **MPI** — distributed-memory parallelism

---

## Platform Support

OSLO has been tested on the following configurations:

| OS             | Compiler       | Status |
| -------------- | -------------- | ------ |
| macOS          | GNU (gfortran) | ✅      |
| Ubuntu Linux   | GNU            | ✅      |
| Ubuntu Linux   | Intel          | ✅      |
| OpenSUSE Linux | GNU            | ✅      |
| OpenSUSE Linux | Intel          | ✅      |

The toolbox is expected to work on most Unix-like systems with a compatible Fortran compiler.

---

## Installation

### Quick Start

```bash
git clone https://github.com/MarcoGrossi92/OSLO.git
cd OSLO
./install.sh build --compiler=gnu
```

### Build Variants

```bash
./install.sh [GLOBAL_OPTIONS] build [OPTIONS]
```

**Global options**

* `-v, --verbose` — verbose output

**Build options**

* `--compiler=gnu|intel`
* `--use-openmp`
* `--use-mpi`
* `--use-sundials` (enabled by default)

**Examples**

```bash
# GNU compiler + OpenMP
./install.sh build --compiler=gnu --use-openmp

# Intel compiler + MPI
./install.sh build --compiler=intel --use-mpi
```

---

## Manual CMake Build

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DUSE_OPENMP=ON \
      -DUSE_MPI=OFF \
      -DUSE_SUNDIALS=ON \
      ..
make
make install
```

`-DOSLO_SUNDIALS_DIR=<folder>` installs SUNDIALS in `<folder>` instead of `sundials-install`. An installation is used
as long as its stamp (`oslo-stamp.txt`: the SUNDIALS commit, the compilers and the flags) matches the build, and is
rebuilt otherwise; the SUNDIALS sources (the submodule `lib/sundials`) are cloned only to build it.

---

## Testing

A full build compiles the test programs in `bin/test/`; run them with `ctest`:

```bash
ctest --test-dir build --output-on-failure
```

Each program checks its results and stops with a nonzero exit status when one is wrong:

| Test                  | What it checks                                                                                       |
| --------------------- | ---------------------------------------------------------------------------------------------------- |
| `OSLO.robertson`      | Every stiff solver on Robertson's problem, against reference values, up to t = 4e10; also with the analytical Jacobian (IJAC = 1) |
| `OSLO.exact_solutions`| Every solver on problems with a known solution: Prothero-Robinson (stiff) and the harmonic oscillator |
| `OSLO.options`        | `solout`, `hmax`, interrupting from `solout`, and the failure status when the step limit is reached   |
| `OSLO.robertsonOMP`   | (with OpenMP) every solver from a parallel region on 4 threads: each result equal to the serial one   |

GitHub Actions (`.github/workflows/tests.yml`) builds OSLO and runs the suite on every pull request into `main`, in
three configurations: gfortran with OpenMP and SUNDIALS, gfortran without either, and ifx with OpenMP and SUNDIALS.

The pre-push hook in `.githooks/` runs the same suite on the build in `build/` and stops the push when a test fails.
Activate it once per clone:

```bash
git config core.hooksPath .githooks
```

---

## Usage

OSLO builds a Fortran library that can be linked into your own applications.

### Minimal Example

```fortran
use OSLO
implicit none

integer, parameter :: neq = 3
real(8) :: y(neq), t, tout
real(8) :: rtol(neq), atol(neq)
integer :: err

call setup_odesolver(N=neq, solver='dvodef90', RT=rtol, AT=atol)
call run_odesolver(neq, t, tout, y, Fgeneral, err)
```

### Linking Against OSLO

```bash
gfortran -c your_program.f90 -Ibuild/modules
gfortran your_program.o -Lbuild/lib -lOSLO -o your_program
```

---

## Third-Party Software

### Bundled

* **Hairer & Wanner Solvers**  
  © University of Geneva — BSD-style license  
  [https://www.unige.ch/~hairer/software.html](https://www.unige.ch/~hairer/software.html)

* **FATODE**  
  Fast and Accurate ODE solver with AD — GNUv1.3 license  
  [https://people.cs.vt.edu/~asandu/Software/FATODE/](https://people.cs.vt.edu/~asandu/Software/FATODE/)

### Optional

* **SUNDIALS**  
  © Lawrence Livermore National Laboratory — BSD license  
  [https://computing.llnl.gov/projects/sundials](https://computing.llnl.gov/projects/sundials)

---

## Project Structure

```text
OSLO/
├── src/            # Core library and solvers
├── lib/            # Third-party solver sources
├── bin/            # Executables and tests
├── cmake/          # CMake modules
├── scripts/        # Helper scripts
├── install.sh      # Build helper
└── CMakeLists.txt
```

---

## Contributing

Contributions are welcome — bug fixes, solver additions, documentation, or performance improvements.

Typical workflow:

1. Fork the repo
2. Create a feature branch
3. Commit changes
4. Open a pull request

---

## Support

Please use the GitHub issue tracker for:

* bug reports
* feature requests
* solver questions

👉 [https://github.com/MarcoGrossi92/OSLO/issues](https://github.com/MarcoGrossi92/OSLO/issues)

---
