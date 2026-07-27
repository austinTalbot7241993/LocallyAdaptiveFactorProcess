#!/usr/bin/env python3
"""
Multivariate Locally Adaptive Factor Process (LAFP) Demonstration.

This script demonstrates fitting a 3-variate continuous time series with
time-varying dynamic covariance structures using the `MultivariateLAFPFactorModel`
(Durante, Scarpa, & Dunson, JMLR 2014).
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from lafp import MultivariateLAFPFactorModel

def main():
    print("==========================================================================")
    print(" Multivariate Locally Adaptive Factor Process (Durante et al. 2014) Demo")
    print("==========================================================================")

    # 1. Synthesize 3-variate time series data with regime shift in volatility & correlation
    np.random.seed(42)
    n_t = 100
    tobs = np.linspace(0.0, 10.0, n_t)

    # True time-varying signals
    mu1 = np.sin(tobs)
    mu2 = np.cos(tobs * 0.8)
    mu3 = 0.5 * np.abs(np.sin(tobs))


    # Time-varying volatility: sharp regime shift around t=5.0
    vol = np.where(tobs > 5.0, 0.4, 0.05)

    y1 = mu1 + vol * np.random.randn(n_t)
    y2 = mu2 + vol * np.random.randn(n_t)
    y3 = mu3 + vol * np.random.randn(n_t)
    y = np.column_stack([y1, y2, y3])

    print(f"Data generated: {n_t} time points, {y.shape[1]} response variables.")
    print("Fitting Multivariate LAFP model (K=2 factors, L=2 dictionary dimension)...")

    # 2. Fit Multivariate LAFP model
    model = MultivariateLAFPFactorModel(
        n_factors=2,
        n_dict=2,
        n_iter=500,  # 500 iterations thinned by 5 => 100 posterior draws
    )
    model.fit(y, tobs)

    print("Model fitting completed successfully!")
    print(f"  Posterior mean draws shape (mu_)    : {model.mu_.shape}")
    print(f"  Posterior cov draws shape (Sigma_) : {model.Sigma_.shape}")

    # Compute posterior mean across MCMC draws
    post_mu_mean = np.mean(model.mu_, axis=0)       # (n_t, 3)
    post_cov_mean = np.mean(model.Sigma_, axis=0)   # (n_t, 3, 3)

    # Extract time-varying variances
    var_y1 = post_cov_mean[:, 0, 0]
    var_y2 = post_cov_mean[:, 1, 1]
    var_y3 = post_cov_mean[:, 2, 2]

    # Extract time-varying covariance between y1 and y2
    cov_y1_y2 = post_cov_mean[:, 0, 1]

    # 3. Generate diagnostic plot
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    # Plot 1: Observed data vs Posterior Mean Trajectories
    axes[0].plot(tobs, y1, 'r.', alpha=0.3, label='$y_1$ (Observed)')
    axes[0].plot(tobs, y2, 'g.', alpha=0.3, label='$y_2$ (Observed)')
    axes[0].plot(tobs, y3, 'b.', alpha=0.3, label='$y_3$ (Observed)')
    axes[0].plot(tobs, post_mu_mean[:, 0], 'r-', linewidth=2, label=r'$\hat{\mu}_1(t)$')
    axes[0].plot(tobs, post_mu_mean[:, 1], 'g-', linewidth=2, label=r'$\hat{\mu}_2(t)$')
    axes[0].plot(tobs, post_mu_mean[:, 2], 'b-', linewidth=2, label=r'$\hat{\mu}_3(t)$')
    axes[0].set_title("Multivariate Time Series Means (Durante et al. 2014)")
    axes[0].set_ylabel("Response")
    axes[0].legend(loc="upper right", ncol=2)
    axes[0].grid(True, linestyle="--", alpha=0.5)

    # Plot 2: Time-Varying Variances (Volatilities)
    axes[1].plot(tobs, var_y1, 'r-', linewidth=2, label=r'Var$(y_1(t))$')
    axes[1].plot(tobs, var_y2, 'g-', linewidth=2, label=r'Var$(y_2(t))$')
    axes[1].plot(tobs, var_y3, 'b-', linewidth=2, label=r'Var$(y_3(t))$')
    axes[1].axvline(x=5.0, color='gray', linestyle=':', label='True Volatility Shift (t=5.0)')
    axes[1].set_title("Estimated Time-Varying Marginal Variances (Locally Adaptive Smoothness)")
    axes[1].set_ylabel("Variance")
    axes[1].legend(loc="upper left")
    axes[1].grid(True, linestyle="--", alpha=0.5)

    # Plot 3: Time-Varying Covariance y1-y2
    axes[2].plot(tobs, cov_y1_y2, 'purple', linewidth=2, label=r'Cov$(y_1(t), y_2(t))$')
    axes[2].set_title("Estimated Dynamic Co-Volatility / Covariance Structure")
    axes[2].set_xlabel("Time (t)")
    axes[2].set_ylabel("Covariance")
    axes[2].legend(loc="upper left")
    axes[2].grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()

    out_plot = os.path.join(os.path.dirname(__file__), "lafp_multivariate_demo_plot.png")
    plt.savefig(out_plot, dpi=150)
    print(f"Diagnostic plot saved to: {out_plot}")

if __name__ == "__main__":
    main()
