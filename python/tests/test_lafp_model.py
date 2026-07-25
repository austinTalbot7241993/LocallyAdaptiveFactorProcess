import pytest
import numpy as np
from lafp import LAFPFactorModel

def test_lafp_model_init():
    model = LAFPFactorModel(n_iter=100, sig_u=500.0, sig_a=2.0)
    assert model.n_iter == 100
    assert model.sig_u == 500.0
    assert model.sig_a == 2.0
    assert not model.is_fitted_
    assert model.theta_ is None
    assert model.sig_ is None

def test_lafp_model_fit():
    # Generate synthetic time series
    np.random.seed(42)
    tobs = np.linspace(0.0, 5.0, 30)
    y = np.sin(tobs) + 0.1 * np.random.randn(30)

    model = LAFPFactorModel(n_iter=50)
    res = model.fit(y, tobs)

    assert res is model
    assert model.is_fitted_
    assert model.theta_ is not None
    assert model.sig_ is not None

    # Check output shapes
    # theta_: (Niter, Nt, 3) = (50, 30, 3)
    assert model.theta_.shape == (50, 30, 3)
    # sig_: (Niter, 3) = (50, 3)
    assert model.sig_.shape == (50, 3)

    # Verify elements are valid finite numbers
    assert np.all(np.isfinite(model.theta_))
    assert np.all(np.isfinite(model.sig_))

def test_lafp_model_dimension_mismatch():
    tobs = np.linspace(0.0, 5.0, 30)
    y = np.sin(tobs)[:20]  # Length 20 vs 30

    model = LAFPFactorModel(n_iter=50)
    with pytest.raises(ValueError, match="Length mismatch"):
        model.fit(y, tobs)
