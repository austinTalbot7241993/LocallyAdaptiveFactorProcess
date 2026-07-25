# Locally Adaptive Factor Process (LAFP)

![CI Status](https://github.com/austinTalbot7241993/LocallyAdaptiveFactorProcess/workflows/CI/badge.svg)

High-performance C implementation of the model described in *"Locally adaptive factor processes for multivariate time series"* ([Durante & Dunson, JMLR 2014](https://www.jmlr.org/papers/volume15/durante14a/durante14a.pdf)).

This repository provides state-space filtering, Fast State Smoothing (FSS), Durbin-Koopman simulation smoothing, and Markov chain Monte Carlo (MCMC) posterior sampling routines for non-Gaussian process factor models.

---

## Performance & Design

Higher-level dynamic languages (Python, R, Julia) incur significant garbage collection and allocation overhead when constructing matrices during iterative MCMC sampling. 

This library avoids allocation overhead by pre-allocating state containers at initialization (`NGPmcmc_construct`). During the sampling iterations:
- **Zero dynamic memory allocations/deallocations** occur inside the MCMC loop.
- Matrix workspaces and GSL buffers are reused in-place across iterations.
- Explicit memory deallocation is executed upon object destruction (`NGPmcmc_free`).

---

## Directory Structure

```
LocallyAdaptiveFactorProcess/
├── include/lafp/       # Public C headers (NGPmcmc.h, SSsimulate2.h, KalmanFilter2.h, etc.)
├── src/                # Library implementation (.c source files)
├── examples/           # Command-line runners (run_ngpmcmc.c)
├── tests/              # CTest test suite (test_ngpmcmc.c, test_heinkel.c)
├── cmake/              # CMake config templates for find_package(lafp)
└── .github/workflows/  # CI pipeline (Linux, macOS, Sanitizers, Pthreads)
```

---

## Building & Installation

### Prerequisites

| Platform | Dependencies | Installation Command |
|---|---|---|
| **macOS** | CMake, Ninja, GSL | `brew install cmake ninja gsl` |
| **Ubuntu / Debian** | CMake, Ninja, GSL, GCC | `sudo apt install cmake ninja-build libgsl-dev gcc` |

### Build Commands

```bash
# Configure (Release build by default)
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release

# Build the library, examples, and test suite
cmake --build build --parallel

# Execute unit and integration tests
cd build && ctest --output-on-failure
```

### Sanitizer / Debug Build (ASan + UBSan)

To test with AddressSanitizer and UndefinedBehaviorSanitizer memory leak detection enabled:

```bash
cmake -S . -B build_asan -G Ninja -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_C_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
cmake --build build_asan
cd build_asan && ctest --output-on-failure
```

---

## CMake Integration (`find_package`)

To use `lafp` in an external CMake project:

```cmake
find_package(lafp REQUIRED)

add_executable(my_app main.c)
target_link_libraries(my_app PRIVATE lafp::lafp_shared)
```

---

## License

This software is released under the MIT License.
