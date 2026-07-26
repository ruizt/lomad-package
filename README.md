# lomad

**Lo**cal **M**oving **A**verage **D**ecoupling — an R package for detecting 
periods of trend similarity and decoupling between two mean-nonstationary time 
series.

## Overview

`lomad` implements a method for detecting localized periods of dissimilarity 
between otherwise correlated moving averages. The method:

1. Smooths each series with a moving average filter to isolate trends.
2. Computes rolling window correlations on the smoothed series.
3. Tests pointwise whether observed correlations fall below their expected
   values under a shared-trend null, using a CLT-based test statistic with
   Benjamini–Yekutieli FDR correction.

Noise parameters (AR/ARMA) are estimated via a variogram-based approach that
is robust to trend contamination (Hall and Van Keilegom, 2003). The asymptotic 
variance of the local correlation accounts for autocorrelation in the smoothed 
noise.

### Core workflow

| Function | Role |
|---|---|
| `lomad_fit()` | Fit the null model (estimate noise, compute expected correlations and asymptotic variance) |
| `lomad_test()` | Pointwise test for local decoupling with FDR correction |
| `lomad()` | Convenience wrapper: `lomad_fit()` → `lomad_test()` |

### Supporting functions

| Function | Role |
|---|---|
| `lomad_plot()` | Visualise fit and test results |
| `estimate_trends()` | Extract trends via moving average |
| `estimate_ar1_noise()` | Estimate AR(1) noise parameters via variogram |
| `estimate_arma_noise()` | Estimate ARMA noise parameters via long-AR approximation |
| `estimate_acf_noise()` | Nonparametric noise autocovariance via the variogram (model-free fallback) |
| `sim_trends()` | Generate synthetic trend pairs at controlled L² separation |
| `sim_noise_pair()` | Add calibrated ARMA noise to trend pairs |

## Installation

### From GitHub

```r
# pak (recommended)
pak::pak("ruizt/lomad-package")

# or, with remotes
remotes::install_github("ruizt/lomad-package")
```

### Development version (from source)

Clone the repository and install with `pak`:

```bash
git clone https://github.com/ruizt/lomad-package.git
cd lomad-package
```

```r
pak::local_install()   # installs dependencies, then the package
```

Or build and install from the command line:

```bash
R CMD build .
R CMD INSTALL lomad_*.tar.gz
```

To load the package in-place during development:

```r
devtools::load_all()
```

## Quickstart

```r
library(lomad)

# Simulate paired trends with controlled L² separation
trends <- sim_trends(n = 500, d = 2, method = "smooth", bw = 50, seed = 1)

# Add calibrated AR(1) noise at target SNR
sim <- sim_noise_pair(trends, h = 10, lambda_target = 1.5,
                      ar.coefs = 0.5, seed = 2)

# Fit null model and test for decoupling
fit <- lomad_fit(sim$y1, sim$y2, h = 10, s = 50)
tst <- lomad_test(fit, alpha = 0.05)

# Which time points show significant decoupling?
which(tst$rejected)

# Visualise
lomad_plot(fit, tst)
```

### Convenience wrapper

```r
out <- lomad(sim$y1, sim$y2, h = 10, s = 50, alpha = 0.05)
out$fit   # lomad_fit output
out$test  # lomad_test output
```

## Citation

If you use `lomad`, please cite the methods paper (and the software where the
specific implementation is relevant):

> Ruiz, T. D., Seifert, A. J., Hamilton, E., Mispagel, C. M., Hunt, O. P.,
> Garcia, J., and Bockmon, E. E. (2026). Inference for local trend similarity
> in nonstationary time series via rolling correlation, with application to
> assessing stability in an estuarine system. Manuscript in preparation.

Run `citation("lomad")` for BibTeX entries for both the paper and the software.
Method references for individual functions appear in their help pages (e.g.
`?estimate_arma_noise`, `?lomad_test`).

## For contributors

Clone the repo and open `lomad.Rproj` in RStudio. Load all package functions
with `devtools::load_all()`.

See `AGENTS.md` for source conventions.

All analyses reported in the methods paper — numerical validation, simulation 
studies, and the Morro Bay field analysis — live in the companion repository
[`lomad-analysis`](https://github.com/ruizt/lomad-analysis), which installs
this package from GitHub.

## Project structure

```
lomad-package/
├── R/                     # Package source (themed: sim_*, estimate_*, lomad_*, utils-*)
├── tests/testthat/        # Unit tests (themed by module)
├── man/                   # Auto-generated documentation
├── inst/CITATION
├── DESCRIPTION
├── NAMESPACE
└── AGENTS.md              # Source conventions for developers and AI agents
```
