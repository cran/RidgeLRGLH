### R code from vignette source 'RidgeLRGLH.Rnw'

###################################################
### code chunk number 1: example
###################################################
library(RidgeLRGLH)
set.seed(8)
Y <- matrix(rnorm(20 * 80), 20)
X <- rbind(1, matrix(rnorm(10 * 80), 10))
C <- diag(11)[, -1, drop = FALSE]
fit <- ridge_glht(Y, X, C, selection = "minimax",
                 control = list(size = 60, grid_size = 80))
fit
head(fit$tuning)


###################################################
### code chunk number 2: diagnostics
###################################################
fit$spectrum$loss
sum(fit$spectrum$weights)
fit$calibration


###################################################
### code chunk number 3: oracle
###################################################
ridge_glht(Y, X, C, lambda = 1,
           spectrum = list(sigma = 1, weights = 1))
