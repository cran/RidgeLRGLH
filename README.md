# RidgeLRGLH

Ridge-regularized largest root tests for high-dimensional general linear
hypotheses, implementing Haoran Li's article accepted by the *Annals of
Statistics*. The public preprint is available at
<https://arxiv.org/abs/2504.15510>.

## Install

From the directory containing the source archive:

```r
install.packages(c("lpSolve", "RMTstat"))
install.packages("RidgeLRGLH_0.1.0.tar.gz", repos = NULL, type = "source")
library(RidgeLRGLH)
```

## Multivariate regression

Rows of `Y` are response variables; columns are independent observations.
Rows of `X` are predictors. Include the intercept explicitly. The hypothesis
is `B %*% C = 0` in `Y = B %*% X + E`.

```r
set.seed(8)
Y <- matrix(rnorm(80 * 150), nrow = 80)
X <- rbind(intercept = 1, matrix(rnorm(20 * 150), nrow = 20))
C <- diag(21)[, -1, drop = FALSE]

fixed <- ridge_glht(Y, X, C)  # ridge defaults to mean residual eigenvalue
adaptive <- ridge_glht(Y, X, C, selection = "minimax")
adaptive
adaptive$tuning
adaptive$spectrum[c("sigma", "weights", "loss")]
```

Other adaptive choices are `"bayes_identity"` and `"bayes_covariance"`.
Set `lambda` for a fixed ridge or `lambda_grid` for adaptive selection.
Do not select the smallest p-value across methods and interpret it as a
single calibrated test. The reported p-values use an asymptotic
Tracy-Widom approximation under the article's assumptions.

## Precomputed covariance matrices

```r
set.seed(17)
Z1 <- matrix(rnorm(40 * 30), 40)
Z2 <- matrix(rnorm(40 * 80), 40)
W1 <- tcrossprod(Z1) / 30
W2 <- tcrossprod(Z2) / 80
ridge_root(W1, W2, n1 = 30, n2 = 80, lambda = 1)
```

The matrices must be normalized by their respective degrees of freedom.
`ridge_root()` is a double Wishart test, not a general equality-of-covariance
matrices test for arbitrary pairs of sample covariance matrices.

## Numerical controls

Defaults use 500 transform points and 500 population-spectrum masses with
first- and second-order spectral constraints. Use `control = list(size = 100,
grid_size = 150)` for exploratory computations. Compare results with larger
grids before drawing conclusions. Inspect `fit$spectrum$loss`; the package
does not assign a universal statistical cutoff to this approximation loss.

Supply `spectrum = list(sigma = ..., weights = ...)` to use a known population
spectrum, for example in oracle experiments. `ridge_calibration()` exposes
centering, scaling, the regular edge and transform derivatives independently.
Read `vignette("RidgeLRGLH", package = "RidgeLRGLH")` for a worked introduction.

## Development checks

```sh
R CMD build RidgeLRGLH
R CMD check --as-cran RidgeLRGLH_0.1.0.tar.gz
```

Tests cover analytic limits, an independent ODE implementation, GLHT projection
identities, singular covariance matrices, adaptive criteria and invalid inputs.
The original manuscript and simulation files remain outside this package.

## Finite-sample validation

Seeded Gaussian identity-covariance checks found 5 rejections in 100 null
replicates for default minimax selection at `p = 200, n1 = 100, n2 = 400`.
At the smaller setting `p = 60, n1 = 40, n2 = 120`, it rejected 26 of 200
null replicates (13%) at nominal 5%, using the full default grids. Increasing
the spectral grid alone did not resolve that finite-sample inflation.
These limited experiments do not establish accuracy across the model class.
Use particular care with smaller dimensions and validate calibration for
your intended setting. The package implements the published asymptotic
method without an additional finite-sample correction.
