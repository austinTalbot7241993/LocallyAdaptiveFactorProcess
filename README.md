# Locally Adaptive Factor Process (LAFP)

![CI Status](https://github.com/austinTalbot7241993/LocallyAdaptiveFactorProcess/workflows/CI/badge.svg)

High-performance C & Python implementation of locally adaptive dynamic factor processes and non-Gaussian process state-space models.

This library implements two primary non-parametric state-space models:
1. **Multivariate LAFP Model**: Dynamic factor processes and dynamic covariance matrix estimation for $p$-dimensional time series ([Durante, Scarpa, & Dunson, JMLR 2014](https://www.jmlr.org/papers/volume15/durante14a/durante14a.pdf)).
2. **Univariate nGP Model**: Continuous time non-Gaussian process regression via nested SDEs for 1D time series ([Zhu & Dunson, JASA 2013](https://arxiv.org/abs/1201.4403)).

---

## Performance & Design

Higher-level dynamic languages (Python, R, Julia) incur significant garbage collection and allocation overhead when constructing state matrices during iterative MCMC sampling. 

This library avoids allocation overhead by pre-allocating state containers at initialization (`NGPlaf2_construct` / `NGPmcmc_construct`). During sampling iterations:
- **Zero dynamic memory allocations/deallocations** occur inside the MCMC loop.
- Matrix workspaces, GSL buffers, and simulation smoothers are reused in-place across iterations.
- Explicit memory deallocation is executed upon object destruction (`NGPlaf2_free` / `NGPmcmc_free`).

---

## Python Package (`lafp`)

The repository includes a native Python package featuring scikit-learn compatible estimator APIs (`MultivariateLAFPFactorModel` and `LAFPFactorModel`).

### Installation

```bash
# Build the C shared library first:
cmake -S . -B build && cmake --build build

# Install the Python package in editable mode:
pip install -e .
```

### Python Quickstart

#### 1. Multivariate LAFP ($p$-dimensional Dynamic Covariance)
```python
import numpy as np
from lafp import MultivariateLAFPFactorModel

# Prepare time vector and 3-variate observation matrix (N_t x 3)
tobs = np.linspace(0.0, 10.0, 100)
y1 = np.sin(tobs) + 0.1 * np.random.randn(100)
y2 = np.cos(tobs) + 0.1 * np.random.randn(100)
y3 = np.sin(tobs) + np.cos(tobs) + 0.1 * np.random.randn(100)
y = np.column_stack([y1, y2, y3])

# Fit multivariate LAFP model (K=2 factors, L=2 dictionary dimension)
model = MultivariateLAFPFactorModel(n_factors=2, n_dict=2, n_iter=1000)
model.fit(y, tobs)

# Access posterior draws for dynamic means and covariances
print("Posterior mean draws shape (mu_)   :", model.mu_.shape)     # (200, 100, 3)
print("Posterior cov draws shape (Sigma_) :", model.Sigma_.shape)  # (200, 100, 3, 3)
```

#### 2. Univariate nGP (1D Time Series)
```python
from lafp import LAFPFactorModel

model = LAFPFactorModel(n_iter=2000, sig_u=1000.0, sig_a=5.0, sig_eps=2.0)
model.fit(y1, tobs)

print("Posterior state draws shape :", model.theta_.shape)  # (2000, 100, 3)
print("Posterior scale draws shape :", model.sig_.shape)    # (2000, 3)
```

---

## Directory Structure

```
LocallyAdaptiveFactorProcess/
├── include/lafp/       # Public C headers (NGPlaf2.h, NGPmcmc.h, SSsimulate2.h, KalmanFilter2.h)
├── src/                # Library implementation (.c source files & lafp-fit CLI)
├── python/             # Python package (lafp/ model.py, _bindings.py, tests/)
├── demos/              # Python & shell demonstration scripts (demo_multivariate_lafp.py, run_multivariate_cli_demo.sh)
├── examples/           # Command-line runners (run_ngpmcmc.c)
├── tests/              # CTest test suite (test_ngplaf2, test_ngpmcmc, test_kalman, test_cli)
├── cmake/              # CMake config templates for find_package(lafp)
└── .github/workflows/  # CI pipeline (Linux, macOS, Sanitizers, Pthreads)
```

---

## Command-Line Interface (`lafp-fit`)

The repository includes a standalone CLI executable (`lafp-fit`) for running MCMC model fitting directly from tabular text files without writing C code.

### Multivariate CLI Example

```bash
./build/lafp-fit -y demos/y_multivariate_demo.txt -t demos/tobs_demo.txt -k 2 -l 2 -n 1000 -o my_mv_experiment
```

Output generated:
- `my_mv_experiment_Mu.txt`: Posterior mean trajectories ($N_{samps} \times N_t \times p$)
- `my_mv_experiment_Sigma0.txt`: Posterior dynamic covariance matrices ($N_{samps} \times N_t \times p \times p$)
- `my_mv_experiment_yhat.txt`: Posterior fitted observation draws
- `my_mv_experiment_Ksi.txt`: Posterior dictionary process draws
- `my_mv_experiment_Psi.txt`: Posterior factor mean draws

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

---

## References

1. Daniele Durante, Bruno Scarpa, and David B. Dunson. (2014). *"Locally Adaptive Factor Processes for Multivariate Time Series"*. **Journal of Machine Learning Research**, 15(44), 1493–1522.
2. Bin Zhu and David B. Dunson. (2013). *"Locally Adaptive Bayes Nonparametric Regression via Nested Gaussian Processes"*. **Journal of the American Statistical Association**, 108(504), 1445–1456.

---

## License

This software is released under the MIT License.
