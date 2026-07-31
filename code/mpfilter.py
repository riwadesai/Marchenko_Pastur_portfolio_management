"""Marchenko-Pastur covariance cleaning.

Canonical implementation used by all phases (full-sample analysis and the
rolling-window backtest). The filtering rule, applied to a correlation matrix C
estimated from T observations of N stocks:

1. Eigendecompose C = V diag(lam) V^T.
2. Determine the noise threshold lambda_+ = sigma2 * (1 + sqrt(N/T))^2.
3. Replace every eigenvalue below lambda_+ by the mean of all such eigenvalues
   (this preserves the trace exactly: variance is reattributed, never created
   or destroyed). Eigenvectors are untouched.
4. Renormalize the reconstructed matrix to unit diagonal,
   C_clean = D^{-1/2} C_f D^{-1/2}, D = diag(C_f), so that rescaling by sample
   volatilities returns each stock's own variance to the diagonal.

Handles the rank-deficient case T < N (q > 1) transparently: the N - T zero
eigenvalues fall below lambda_+ and are lifted to the common noise average,
which is precisely what restores invertibility.
"""
import numpy as np


def mp_lambda_plus(N, T, sigma2):
    return sigma2 * (1.0 + np.sqrt(N / T)) ** 2


def clean_correlation(C, T, sigma2="market"):
    """MP-filter a correlation matrix. Returns (C_clean, info).

    sigma2: "market" -> 1 - lam_max/N (per-window, no fitting required),
            "naive"  -> 1.0,
            float    -> explicit value (e.g. a KS-fitted sigma2).
    """
    C = np.asarray(C, float)
    N = C.shape[0]
    lam, V = np.linalg.eigh(C)                       # ascending
    if sigma2 == "market":
        s2 = 1.0 - lam[-1] / N
    elif sigma2 == "naive":
        s2 = 1.0
    else:
        s2 = float(sigma2)
    lp = mp_lambda_plus(N, T, s2)

    noise = lam < lp
    lam_clean = lam.copy()
    if noise.any():
        lam_clean[noise] = lam[noise].mean()
    C_f = (V * lam_clean) @ V.T

    d = np.sqrt(np.diag(C_f))
    C_clean = C_f / np.outer(d, d)
    np.fill_diagonal(C_clean, 1.0)

    info = {
        "sigma2": float(s2),
        "lambda_plus": float(lp),
        "n_noise": int(noise.sum()),
        "n_signal": int(N - noise.sum()),
        "noise_mean": float(lam[noise].mean()) if noise.any() else None,
        "max_diag_dev_pre_renorm": float(np.abs(d**2 - 1).max()),
        "trace_after_filter": float(lam_clean.sum()),
    }
    return C_clean, info


def clean_covariance(returns, sigma2="market"):
    """T x N returns window -> (Sigma_raw, Sigma_clean, info).

    Standardizes columns (ddof=0), cleans the correlation matrix, and rescales
    both raw and cleaned correlation by the window's sample volatilities.
    """
    X = np.asarray(returns, float)
    T = X.shape[0]
    sd = X.std(0, ddof=0)
    Xs = (X - X.mean(0)) / sd
    C = (Xs.T @ Xs) / T
    C_clean, info = clean_correlation(C, T, sigma2)
    scale = np.outer(sd, sd)
    return C * scale, C_clean * scale, info
