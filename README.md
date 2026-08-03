# From Wigner's Surmise to the Marchenko–Pastur Law

**An Account of Random Matrix Theory, with an Application to NSE Equity Portfolios**

📄 **[Read the paper (PDF)](main.pdf)**

This is an expository paper: a self-contained walk through classical random
matrix theory — what a random matrix is, why the Gaussian ensembles (GOE, GUE)
sit exactly at the intersection of *independent entries* and *rotational
invariance* (proved by a two-by-two trace computation), Wigner's surmise for
eigenvalue spacings (derived in full via a Jacobian change of variables), the
Wishart ensemble, and a complete derivation of the Marchenko–Pastur law using
the Stieltjes transform, the Sherman–Morrison update, Sylvester's determinant
identity, and concentration of quadratic forms — followed by an application to
real data.

The application tests one concrete claim on daily returns of **548 National
Stock Exchange of India equities over 1,236 trading days** (July 2021 – June
2026): *the eigenvalues of the sample correlation matrix that fall inside the
Marchenko–Pastur band are genuinely noise.* The evidence comes in three layers:

1. **Distributed like noise** — the eigenvalue bulk matches the MP density.
2. **Interact like noise** — the level spacings match Wigner's surmise and
   decisively reject the Poisson alternative.
3. **Price like noise** — in an out-of-sample portfolio backtest, erasing the
   noise band produces minimum-variance portfolios with **7.8%** annualized
   volatility, against **9.3%** for Ledoit–Wolf shrinkage, **17.7%** for equal
   weights, and **19.4%** for the raw covariance matrix.

## Repository layout

```
main.pdf        the paper
main.tex        LaTeX source
figures/        the five figures used by the paper
code/           the analysis behind the application section
```

The `code/` folder contains the seven Jupyter notebooks (with their outputs
preserved, so they are readable on GitHub without running anything) plus the
canonical filter implementation:

| File | What it does |
|---|---|
| `phase0_data_acquisition.ipynb` | Downloads and cleans the NSE universe (incl. corporate-action verification) |
| `phase1_noise_problem.ipynb` | Correlation matrix, eigenvalue spectrum, the noise problem |
| `phase2_mp_validation.ipynb` | Marchenko–Pastur density fit with Monte Carlo–calibrated KS tests |
| `phase3_goe_spacing.ipynb` | Wigner-surmise spacing statistics of the noise band |
| `phase4_mp_filtering.ipynb` | The eigenvalue-flattening filter and its sanity checks |
| `phase5_randomized_svd.ipynb` | Randomized SVD implementation and benchmarks |
| `phase67_backtest.ipynb` | The out-of-sample portfolio backtest |
| `mpfilter.py` | The MP cleaning filter (imported by the notebooks) |

tly from the frozen snapshot used in the
paper; the notebooks' saved outputs record the exact numbers reported there.
