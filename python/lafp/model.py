import numpy as np
from ._bindings import LAFPBindings

class LAFPFactorModel:
    """
    Locally Adaptive Factor Process (LAFP) Non-Gaussian Process MCMC Sampler.

    Parameters
    ----------
    n_iter : int, default=2000
        Total number of MCMC sampling iterations to execute.
    sig_u : float, default=1000.0
        Initial process noise scale for state components.
    sig_a : float, default=5.0
        Initial process noise scale for derivative components.
    sig_eps : float, default=2.0
        Initial observation noise standard deviation.
    sig_mu : float, default=4.0
        Prior standard deviation for mean state components.
    sig_alpha : float, default=4.0
        Prior standard deviation for derivative state components.
    prior_a : float, default=1.0
        Inverse-Gamma prior shape parameter for noise variance.
    prior_b : float, default=1.0
        Inverse-Gamma prior scale parameter for noise variance.

    Attributes
    ----------
    theta_ : ndarray of shape (n_iter, n_t, 3)
        Posterior state draws sampled during MCMC iterations.
    sig_ : ndarray of shape (n_iter, 3)
        Posterior process scale parameter draws sampled during MCMC iterations.
    is_fitted_ : bool
        Boolean flag indicating whether `.fit()` has been executed.

    Examples
    --------
    >>> import numpy as np
    >>> from lafp import LAFPFactorModel
    >>> tobs = np.linspace(0.0, 10.0, 50)
    >>> y = np.sin(tobs).reshape(-1, 1)
    >>> model = LAFPFactorModel(n_iter=100)
    >>> model.fit(y, tobs)
    >>> model.theta_.shape
    (100, 50, 3)
    """

    def __init__(
        self,
        n_iter: int = 2000,
        sig_u: float = 1000.0,
        sig_a: float = 5.0,
        sig_eps: float = 2.0,
        sig_mu: float = 4.0,
        sig_alpha: float = 4.0,
        prior_a: float = 1.0,
        prior_b: float = 1.0,
    ):
        self.n_iter = n_iter
        self.sig_u = sig_u
        self.sig_a = sig_a
        self.sig_eps = sig_eps
        self.sig_mu = sig_mu
        self.sig_alpha = sig_alpha
        self.prior_a = prior_a
        self.prior_b = prior_b

        self.theta_ = None
        self.sig_ = None
        self.is_fitted_ = False

    def fit(self, y: np.ndarray, tobs: np.ndarray):
        """
        Fit the LAFP MCMC model on time series observation data.

        Parameters
        ----------
        y : array-like of shape (n_t,) or (n_t, n_p)
            Observation values across time points.
        tobs : array-like of shape (n_t,)
            Observation timestamps.

        Returns
        -------
        self : object
            Returns self fitted instance.
        """
        y_arr = np.asarray(y, dtype=np.float64)
        if y_arr.ndim == 1:
            y_arr = y_arr.reshape(-1, 1)

        tobs_arr = np.asarray(tobs, dtype=np.float64).ravel()
        n_t = len(tobs_arr)

        if len(y_arr) != n_t:
            raise ValueError(f"Length mismatch: y has {len(y_arr)} rows, but tobs has {n_t} elements.")

        b = LAFPBindings()
        lib = b.lib

        # Allocate GSL structures
        gsl_y = b.numpy_to_gsl_matrix(y_arr)
        gsl_tobs = b.numpy_to_gsl_vector(tobs_arr)
        gsl_sigu = b.numpy_to_gsl_vector([self.sig_u])
        gsl_siga = b.numpy_to_gsl_vector([self.sig_a])

        # Allocate Output structures
        marray_th = lib.marray3d_alloc(self.n_iter, n_t, 3)
        gsl_sig = lib.gsl_matrix_alloc(self.n_iter, 3)

        # Allocate Sampler Instance
        mcmc = lib.NGPmcmc_New()
        if not mcmc:
            raise RuntimeError("Failed to allocate NGPmcmc sampler memory.")

        try:
            # Construct sampler workspace
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

            # Run MCMC iterations
            ret_ops = lib.NGPmcmc_operations(mcmc)
            if ret_ops != 0:
                raise RuntimeError(f"NGPmcmc_operations failed with return code {ret_ops}")

            # Extract results into NumPy arrays
            self.theta_ = b.marray3d_to_numpy(marray_th, self.n_iter, n_t, 3)
            self.sig_ = b.gsl_matrix_to_numpy(gsl_sig, self.n_iter, 3)
            self.is_fitted_ = True

        finally:
            # Clean up C memory handles
            lib.NGPmcmc_free(mcmc)
            lib.gsl_matrix_free(gsl_y)
            lib.gsl_vector_free(gsl_tobs)
            lib.gsl_vector_free(gsl_sigu)
            lib.gsl_vector_free(gsl_siga)
            lib.marray3d_free(marray_th)
            lib.gsl_matrix_free(gsl_sig)

        return self
