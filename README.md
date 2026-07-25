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

## Command-Line Interface (`lafp-fit`)

The repository includes a standalone CLI executable (`lafp-fit`) for running MCMC model fitting directly from tabular text files without writing C code.

### Quickstart Example

After building the project, run `lafp-fit` on sample data:

```bash
# Basic run with auto-detected length:
./build/lafp-fit -y tests/y.txt -t tests/tobs.txt -n 2000 -o my_experiment

# Output printed to terminal:
# [LAFP-FIT] Starting MCMC Fit...
#   Input y      : tests/y.txt
#   Input tobs   : tests/tobs.txt
#   Time Points  : 1001
#   Iterations   : 2000
#   Out Prefix   : my_experiment
# [LAFP-FIT] Executing MCMC sampling iterations...
# [LAFP-FIT] Exporting output files with prefix 'my_experiment'...
# [LAFP-FIT] MCMC Fit Completed Successfully.
```

The command generates two posterior output text files:
- `my_experiment_Theta.txt`: Posterior state estimates matrix ($N_{iter} \times N_t \times 3$)
- `my_experiment_Sig.txt`: Posterior process scale parameter draws ($N_{iter} \times 3$)

---

### Customizing Priors and Initial Values

Custom noise scales and Inverse-Gamma prior parameters can be passed directly via command-line arguments:

```bash
./build/lafp-fit \
  --input-y tests/y.txt \
  --input-t tests/tobs.txt \
  --niter 5000 \
  --out-prefix run_prior_test \
  --sig-u 500.0 \
  --sig-a 2.5 \
  --sig-eps 1.5 \
  --sig-mu 2.0 \
  --sig-alpha 2.0 \
  --prior-a 2.0 \
  --prior-b 1.0
```

---

### CLI Reference Guide

| Option | Long Flag | Description | Default Value |
|---|---|---|---|
| `-y` | `--input-y <path>` | Path to observation matrix file ($y$) | **Required** |
| `-t` | `--input-t <path>` | Path to time vector file ($t_{obs}$) | **Required** |
| `-n` | `--niter <int>` | Total MCMC sampling iterations | `2000` |
| `-o` | `--out-prefix <str>` | Output filename prefix | `lafp_out` |
| | `--nt <int>` | Number of time points | Auto-detected |
| | `--sig-u <val>` | Initial state noise scale ($\sigma_U$) | `1000.0` |
| | `--sig-a <val>` | Initial derivative noise scale ($\sigma_A$) | `5.0` |
| | `--sig-eps <val>` | Initial observation noise ($\sigma_{\epsilon}$) | `2.0` |
| | `--sig-mu <val>` | Prior scale for mean state | `4.0` |
| | `--sig-alpha <val>` | Prior scale for derivative state | `4.0` |
| `-a` | `--prior-a <val>` | Inverse-Gamma prior shape parameter ($a$) | `1.0` |
| `-b` | `--prior-b <val>` | Inverse-Gamma prior scale parameter ($b$) | `1.0` |
| `-h` | `--help` | Display usage help menu and exit | |
| `-v` | `--version` | Display version information (`1.0.0`) | |

---

## Building & Installation

### Prerequisites

| Platform | Dependencies | Installation Command |
|---|---|---|
| **macOS** | CMake, GSL | `brew install cmake gsl` |
| **Ubuntu / Debian** | CMake, GSL, GCC | `sudo apt install cmake libgsl-dev gcc` |

### Build Commands

```bash
# Configure build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# Build library, lafp-fit CLI, and test suite
cmake --build build --parallel

# Execute unit and integration tests
cd build && ctest --output-on-failure
```

### Sanitizer / Debug Build (ASan + UBSan)

To test with AddressSanitizer and UndefinedBehaviorSanitizer memory leak detection enabled:

```bash
cmake -S . -B build_asan -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_C_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
cmake --build build_asan
cd build_asan && ctest --output-on-failure
```

---

## CMake Integration (`find_package`)

To use `lafp` as a library in an external CMake project:

```cmake
find_package(lafp REQUIRED)

add_executable(my_app main.c)
target_link_libraries(my_app PRIVATE lafp::lafp_shared)
```

---

## License

This software is released under the MIT License.
