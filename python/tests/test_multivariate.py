import numpy as np
import pytest
from lafp import LAFPFactorModel, MultivariateLAFPFactorModel

def test_multivariate_lafp_model_basic():
    """Test basic multivariate LAFP fitting with MultivariateLAFPFactorModel."""
    np.random.seed(42)
    n_t = 50
    n_p = 3
    tobs = np.linspace(0.0, 10.0, n_t)
    
    # Generate 3-variate linearly independent time series
    y = np.column_stack([
        np.sin(tobs) + 0.1 * np.random.randn(n_t),
        np.cos(0.8 * tobs) + 0.1 * np.random.randn(n_t),
        0.5 * np.abs(np.sin(tobs)) + 0.1 * np.random.randn(n_t),
    ])

    model = MultivariateLAFPFactorModel(
        n_factors=2,
        n_dict=2,
        n_iter=500,
    )
    model.fit(y, tobs)

    assert model.is_fitted_ is True
    # 500 iterations thinned by 5 yields 100 saved samples
    n_samps = 100
    assert model.mu_.shape == (n_samps, n_t, n_p)
    assert model.Sigma_.shape == (n_samps, n_t, n_p, n_p)
    assert model.yhat_.shape == (n_samps, n_t, n_p)
    assert model.Ksi_.shape == (n_samps, n_t, 4)  # 2 * 2 = 4
    assert model.Psi_.shape == (n_samps, n_t, 2)

    # Check non-NaN
    assert not np.isnan(model.mu_).any()
    assert not np.isnan(model.Sigma_).any()
    assert not np.isnan(model.yhat_).any()


def test_lafp_factor_model_auto_dispatch():
    """Test LAFPFactorModel auto-dispatching 2D data to multivariate LAFP engine."""
    np.random.seed(123)
    n_t = 50
    tobs = np.linspace(0.0, 10.0, n_t)
    y_2d = np.column_stack([
        np.sin(tobs) + 0.1 * np.random.randn(n_t),
        np.cos(0.8 * tobs) + 0.1 * np.random.randn(n_t),
    ])

    model = LAFPFactorModel(n_iter=500, n_factors=2, n_dict=2)
    model.fit(y_2d, tobs)

    assert model.is_fitted_ is True
    assert model.mu_ is not None
    assert model.Sigma_ is not None
    assert model.mu_.shape == (100, n_t, 2)
    assert model.Sigma_.shape == (100, n_t, 2, 2)
    assert not np.isnan(model.mu_).any()
    assert not np.isnan(model.Sigma_).any()


def test_lafp_univariate_dispatch():
    """Test LAFPFactorModel univariate dispatch for 1D input."""
    np.random.seed(7)
    n_t = 20
    tobs = np.linspace(0.0, 2.0, n_t)
    y_1d = np.sin(tobs)

    model = LAFPFactorModel(n_iter=20)
    model.fit(y_1d, tobs)

    assert model.is_fitted_ is True
    assert model.theta_ is not None
    assert model.sig_ is not None
    assert model.theta_.shape == (20, n_t, 3)
    assert model.sig_.shape == (20, 3)
