import numpy as np
from ._bindings import LAFPBindings

class MultivariateLAFPFactorModel:
    """
    Multivariate Locally Adaptive Factor Process (LAFP) Sampler.

    Implements the dynamic factor and time-varying covariance model from
    Durante, Scarpa, & Dunson (JMLR 2014), "Locally adaptive factor processes for
    multivariate time series".

    Parameters
    ----------
    n_factors : int, default=2
        Number of latent factors (K).
    n_dict : int, default=2
        Number of dictionary process elements (L).
    n_iter : int, default=2000
        Total MCMC sampling iterations to run. Note that posterior output draws
        are thinned by 5, yielding `n_iter // 5` saved samples.
    sig_mu : float, default=4.0
        Prior standard deviation for mean state initial components.
    sig_alpha : float, default=4.0
        Prior standard deviation for derivative initial components.
    prior_eps : tuple of float, default=(1.0, 0.1)
        Inverse-Gamma prior (shape, scale) for observation error variances.
    prior_ksi : tuple of float, default=(1.0, 1.0)
        Inverse-Gamma prior (shape, scale) for dictionary process variance scale.
    prior_a : tuple of float, default=(1.0, 1.0)
        Inverse-Gamma prior (shape, scale) for dictionary mean process scale.
    prior_psi : tuple of float, default=(1.0, 1.0)
        Inverse-Gamma prior (shape, scale) for factor process variance scale.
    prior_b : tuple of float, default=(1.0, 1.0)
        Inverse-Gamma prior (shape, scale) for factor mean process scale.
    prior_aa : tuple of float, default=(3.0, 3.0)
        Gamma prior parameters (a1, a2) for shrinkage factor multipliers.

    Attributes
    ----------
    mu_ : ndarray of shape (n_samples, n_t, n_p)
        Posterior draws of time-varying mean vectors mu(t).
    Sigma_ : ndarray of shape (n_samples, n_t, n_p, n_p)
        Posterior draws of time-varying covariance matrices Sigma(t).
    yhat_ : ndarray of shape (n_samples, n_t, n_p)
        Posterior draws of fitted observation vectors yhat(t).
    Ksi_ : ndarray of shape (n_samples, n_t, n_dict * n_factors)
        Posterior draws of dictionary process matrix Xi(t).
    Psi_ : ndarray of shape (n_samples, n_t, n_factors)
        Posterior draws of factor mean process vector Psi(t).
    is_fitted_ : bool
        True if model has been fitted.

    References
    ----------
    Daniele Durante, Bruno Scarpa, and David B. Dunson. (2014).
    Locally Adaptive Factor Processes for Multivariate Time Series.
    Journal of Machine Learning Research, 15(44), 1493-1522.
    """

    def __init__(
        self,
        n_factors: int = 2,
        n_dict: int = 2,
        n_iter: int = 2000,
        sig_mu: float = 4.0,
        sig_alpha: float = 4.0,
        prior_eps: tuple = (1.0, 0.1),
        prior_ksi: tuple = (1.0, 1.0),
        prior_a: tuple = (1.0, 1.0),
        prior_psi: tuple = (1.0, 1.0),
        prior_b: tuple = (1.0, 1.0),
        prior_aa: tuple = (3.0, 3.0),
    ):
        self.n_factors = n_factors
        self.n_dict = n_dict
        self.n_iter = n_iter
        self.sig_mu = sig_mu
        self.sig_alpha = sig_alpha
        self.prior_eps = prior_eps
        self.prior_ksi = prior_ksi
        self.prior_a = prior_a
        self.prior_psi = prior_psi
        self.prior_b = prior_b
        self.prior_aa = prior_aa

        self.mu_ = None
        self.Sigma_ = None
        self.yhat_ = None
        self.Ksi_ = None
        self.Psi_ = None
        self.is_fitted_ = False

    def fit(self, y: np.ndarray, tobs: np.ndarray):
        """
        Fit multivariate LAFP model to observation data Y (n_t x n_p).
        """
        y_arr = np.asarray(y, dtype=np.float64)
        if y_arr.ndim == 1:
            y_arr = y_arr.reshape(-1, 1)

        n_t, n_p = y_arr.shape
        tobs_arr = np.asarray(tobs, dtype=np.float64).ravel()

        if len(tobs_arr) != n_t:
            raise ValueError(f"Length mismatch: y has {n_t} rows, but tobs has {len(tobs_arr)} elements.")

        b = LAFPBindings()
        lib = b.lib

        n_samps = self.n_iter // 5
        if n_samps < 1:
            raise ValueError(f"n_iter must be at least 5 (got {self.n_iter}).")

        def _to_pair(val):
            if isinstance(val, (int, float)):
                return (float(val), float(val))
            return tuple(val)

        # GSL Inputs
        gsl_y = b.numpy_to_gsl_matrix(y_arr)
        gsl_tobs = b.numpy_to_gsl_vector(tobs_arr)

        gsl_sig_prior = b.numpy_to_gsl_vector([self.sig_mu, self.sig_alpha])
        gsl_eps_prior = b.numpy_to_gsl_vector(_to_pair(self.prior_eps))
        gsl_ksi_prior = b.numpy_to_gsl_vector(_to_pair(self.prior_ksi))
        gsl_a_prior = b.numpy_to_gsl_vector(_to_pair(self.prior_a))
        gsl_psi_prior = b.numpy_to_gsl_vector(_to_pair(self.prior_psi))
        gsl_b_prior = b.numpy_to_gsl_vector(_to_pair(self.prior_b))
        gsl_aa_prior = b.numpy_to_gsl_vector(_to_pair(self.prior_aa))

        # Output marrays
        marray_ksi = lib.marray3d_alloc(n_samps, n_t, self.n_factors * self.n_dict)
        marray_psi = lib.marray3d_alloc(n_samps, n_t, self.n_factors)
        marray_yhat = lib.marray3d_alloc(n_samps, n_t, n_p)
        marray_mu = lib.marray3d_alloc(n_samps, n_t, n_p)
        marray_sigma = lib.marray4d_alloc(n_samps, n_t, n_p, n_p)

        sampler = lib.NGPlaf2_New()
        if not sampler:
            raise RuntimeError("Failed to allocate NGPlaf2 sampler memory.")

        try:
            ret_c = lib.NGPlaf2_construct(
                sampler,
                gsl_tobs,
                gsl_y,
                int(self.n_factors),
                int(self.n_dict),
                int(self.n_iter),
                marray_ksi,
                marray_psi,
                marray_yhat,
                marray_mu,
                marray_sigma,
                gsl_sig_prior,
                gsl_eps_prior,
                gsl_ksi_prior,
                gsl_a_prior,
                gsl_psi_prior,
                gsl_b_prior,
                gsl_aa_prior,
            )
            if ret_c != 0:
                raise RuntimeError(f"NGPlaf2_construct failed with return code {ret_c}")

            ret_ops = lib.NGPlaf2_operations(sampler)
            if ret_ops != 0:
                raise RuntimeError(f"NGPlaf2_operations failed with return code {ret_ops}")

            # Extract arrays
            self.mu_ = b.marray3d_to_numpy(marray_mu, n_samps, n_t, n_p)
            self.yhat_ = b.marray3d_to_numpy(marray_yhat, n_samps, n_t, n_p)
            self.Ksi_ = b.marray3d_to_numpy(marray_ksi, n_samps, n_t, self.n_factors * self.n_dict)
            self.Psi_ = b.marray3d_to_numpy(marray_psi, n_samps, n_t, self.n_factors)
            self.Sigma_ = b.marray4d_to_numpy(marray_sigma, n_samps, n_t, n_p, n_p)
            self.is_fitted_ = True

        finally:
            lib.NGPlaf2_free(sampler)
            lib.gsl_matrix_free(gsl_y)
            lib.gsl_vector_free(gsl_tobs)
            lib.gsl_vector_free(gsl_sig_prior)
            lib.gsl_vector_free(gsl_eps_prior)
            lib.gsl_vector_free(gsl_ksi_prior)
            lib.gsl_vector_free(gsl_a_prior)
            lib.gsl_vector_free(gsl_psi_prior)
            lib.gsl_vector_free(gsl_b_prior)
            lib.gsl_vector_free(gsl_aa_prior)
            lib.marray3d_free(marray_ksi)
            lib.marray3d_free(marray_psi)
            lib.marray3d_free(marray_yhat)
            lib.marray3d_free(marray_mu)
            lib.marray4d_free(marray_sigma)

        return self


class LAFPFactorModel:
    """
    Locally Adaptive Factor Process (LAFP) Estimator.

    Automatically dispatches to:
    - **Multivariate LAFP** (Durante, Scarpa, & Dunson, JMLR 2014) when `n_p > 1` or when `n_factors` is specified.
    - **Univariate nGP** (Zhu & Dunson, JASA 2013) when `n_p == 1` and `n_factors=None`.
    """

    def __init__(
        self,
        n_iter: int = 2000,
        n_factors: int = None,
        n_dict: int = 2,
        sig_u: float = 1000.0,
        sig_a: float = 5.0,
        sig_eps: float = 2.0,
        sig_mu: float = 4.0,
        sig_alpha: float = 4.0,
        prior_a: float = 1.0,
        prior_b: float = 1.0,
    ):
        self.n_iter = n_iter
        self.n_factors = n_factors
        self.n_dict = n_dict
        self.sig_u = sig_u
        self.sig_a = sig_a
        self.sig_eps = sig_eps
        self.sig_mu = sig_mu
        self.sig_alpha = sig_alpha
        self.prior_a = prior_a
        self.prior_b = prior_b

        self.theta_ = None
        self.sig_ = None
        self.mu_ = None
        self.Sigma_ = None
        self.yhat_ = None
        self.is_fitted_ = False
        self._mv_model = None

    def fit(self, y: np.ndarray, tobs: np.ndarray):
        """
        Fit the model on observation matrix `y` and timestamps `tobs`.
        """
        y_arr = np.asarray(y, dtype=np.float64)
        if y_arr.ndim == 1:
            y_arr = y_arr.reshape(-1, 1)

        n_t, n_p = y_arr.shape

        # Use Multivariate model if p > 1 or n_factors is specified
        if n_p > 1 or self.n_factors is not None:
            k = self.n_factors if self.n_factors is not None else min(2, n_p)
            l = self.n_dict if self.n_dict is not None else 2
            self._mv_model = MultivariateLAFPFactorModel(
                n_factors=k,
                n_dict=l,
                n_iter=self.n_iter,
                sig_mu=self.sig_mu,
                sig_alpha=self.sig_alpha,
                prior_a=self.prior_a,
                prior_b=self.prior_b,
            )
            self._mv_model.fit(y_arr, tobs)
            self.mu_ = self._mv_model.mu_
            self.Sigma_ = self._mv_model.Sigma_
            self.yhat_ = self._mv_model.yhat_
            self.Ksi_ = self._mv_model.Ksi_
            self.Psi_ = self._mv_model.Psi_
            self.is_fitted_ = True
            return self

        # Otherwise univariate nGP
        tobs_arr = np.asarray(tobs, dtype=np.float64).ravel()
        if len(y_arr) != len(tobs_arr):
            raise ValueError(f"Length mismatch: y has {len(y_arr)} rows, but tobs has {len(tobs_arr)} elements.")

        b = LAFPBindings()
        lib = b.lib

        gsl_y = b.numpy_to_gsl_matrix(y_arr)
        gsl_tobs = b.numpy_to_gsl_vector(tobs_arr)
        gsl_sigu = b.numpy_to_gsl_vector([self.sig_u])
        gsl_siga = b.numpy_to_gsl_vector([self.sig_a])

        marray_th = lib.marray3d_alloc(self.n_iter, n_t, 3)
        gsl_sig = lib.gsl_matrix_alloc(self.n_iter, 3)

        mcmc = lib.NGPmcmc_New()
        if not mcmc:
            raise RuntimeError("Failed to allocate NGPmcmc sampler memory.")

        try:
            ret_c = lib.NGPmcmc_construct(
                mcmc,
                gsl_y,
                gsl_tobs,
                self.n_iter,
                gsl_sigu,
                gsl_siga,
                float(self.sig_eps),
                float(self.sig_mu),
                float(self.sig_alpha),
                float(self.prior_a),
                float(self.prior_b),
                marray_th,
                gsl_sig,
            )
            if ret_c != 0:
                raise RuntimeError(f"NGPmcmc_construct failed with return code {ret_c}")

            ret_ops = lib.NGPmcmc_operations(mcmc)
            if ret_ops != 0:
                raise RuntimeError(f"NGPmcmc_operations failed with return code {ret_ops}")

            self.theta_ = b.marray3d_to_numpy(marray_th, self.n_iter, n_t, 3)
            self.sig_ = b.gsl_matrix_to_numpy(gsl_sig, self.n_iter, 3)
            self.is_fitted_ = True

        finally:
            lib.NGPmcmc_free(mcmc)
            lib.gsl_matrix_free(gsl_y)
            lib.gsl_vector_free(gsl_tobs)
            lib.gsl_vector_free(gsl_sigu)
            lib.gsl_vector_free(gsl_siga)
            lib.marray3d_free(marray_th)
            lib.gsl_matrix_free(gsl_sig)

        return self
