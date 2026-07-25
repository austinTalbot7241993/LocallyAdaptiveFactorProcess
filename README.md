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
├── src/                # Library implementation (.c source files & lafp-fit CLI)
├── examples/           # Command-line runners (run_ngpmcmc.c)
├── tests/              # CTest test suite (test_ngpmcmc, test_kalman, test_complex, test_cli)
├── cmake/              # CMake config templates for find_package(lafp)
└── .github/workflows/  # CI pipeline (Linux, macOS, Sanitizers, Pthreads)
```

---

## Command-Line Tool (`lafp-fit`)

The repository includes a standalone CLI application (`lafp-fit`) for running MCMC model fitting directly from data files.

```bash
# Basic Usage:
lafp-fit -y <y_file> -t <tobs_file> [options]

# Example: Run 5000 iterations and output files with prefix 'my_results'
lafp-fit -y data/y.txt -t data/tobs.txt -n 5000 -o my_results

# Display help menu and available hyperparameter options:
lafp-fit --help
```

### CLI Options Summary

| Flag | Argument | Description | Default |
|---|---|---|---|
| `-y, --input-y` | `<path>` | Path to observation matrix file | **Required** |
| `-t, --input-t` | `<path>` | Path to time observations file | **Required** |
| `-n, --niter` | `<int>` | Total MCMC sampling iterations | `2000` |
| `-o, --out-prefix` | `<str>` | Prefix for posterior output files | `lafp_out` |
| `--nt` | `<int>` | Number of time points | Auto-detected |
| `--sig-u` | `<float>` | State process noise scale | `1000.0` |
| `--sig-a` | `<float>` | Derivative process noise scale | `5.0` |
| `--sig-eps` | `<float>` | Observation noise std dev | `2.0` |
| `-a, --prior-a` | `<float>` | Inverse-Gamma prior shape parameter | `1.0` |
| `-b, --prior-b` | `<float>` | Inverse-Gamma prior scale parameter | `1.0` |

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

# Build the library, CLI tool, examples, and test suite
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
