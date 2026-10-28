library(RidgeLRGLH)
near <- function(x, y, tol = 1e-7) stopifnot(isTRUE(all.equal(unname(x), unname(y), tolerance = tol, check.attributes = FALSE)))
errors <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), "error"))

# Analytic Marchenko-Pastur limit as denominator degrees of freedom diverge.
for (g in c(0.2, 1, 3)) {
    cal <- ridge_calibration(1, g, 1e-9, 0.7)
    near(cal["center"], (1 + sqrt(g))^2 / 1.7)
    near(cal["scale"], sqrt(g) * (1 + sqrt(g))^(4/3) / 1.7)
    near(cal["beta"]^2 * cal["s_prime"], 1 / g)
}
# Joint covariance/ridge rescaling and repeated spectral atoms.
a <- ridge_calibration(c(1, 2), 0.6, 1.5, 0.5, c(0.4, 0.6))
near(a, ridge_calibration(c(10, 20), 0.6, 1.5, 5, c(0.4, 0.6)))
near(a, ridge_calibration(c(1, 2, 2), 0.6, 1.5, 0.5, c(0.4, 0.3, 0.3)))
# Independent RK4 integration of the ODE, evaluated at the scalar solution.
# This compares both s and the manuscript's derivative formulas.
sigma <- c(0.3, 1, 2)
w <- c(0.2, 0.5, 0.3)
g1 <- 0.7; g2 <- 1.4; lam <- 0.8
cal <- ridge_calibration(sigma, g1, g2, lam, w)
H <- function(h, j) sum(w * (sigma / (lam - sigma * h))^j)
s0 <- uniroot(function(s) H(-1 / (1 + g2 * s), 1) - s, c(0, 10), tol = 1e-12)$root
rhs <- function(x, s) {
    H2 <- H(x - 1 / (1 + g2 * s), 2)
    H2 / (1 - g2 * H2 / (1 + g2 * s)^2)
}
dt <- unname(cal["beta"]) / 2000
s <- s0
for (i in 0:1999) {
    x <- i * dt
    k1 <- rhs(x, s)
    k2 <- rhs(x + dt / 2, s + dt * k1 / 2)
    k3 <- rhs(x + dt / 2, s + dt * k2 / 2)
    k4 <- rhs(x + dt, s + dt * k3)
    s <- s + dt * (k1 + 2*k2 + 2*k3 + k4) / 6
}
near(s, cal["s"], 1e-8)
sp <- rhs(cal["beta"], s)
near(sp, cal["s_prime"], 1e-8)
zeta <- g2 / (1 + g2*s)^2
h <- cal["beta"] - 1/(1 + g2*s)
spp <- (2 * (1 + zeta*sp)^2 * H(h, 3) -
    2 * (1 + g2*s)^(-3) * H(h, 2) * (g2*sp)^2) / (1 - zeta*H(h, 2))
near(spp, cal["s_double_prime"], 1e-8)

# GLHT QR construction versus independently formed projection matrices.
set.seed(110)
p <- 12; n <- 60; m <- 9
Y <- matrix(rnorm(p*n), p)
X <- rbind(1, matrix(rnorm((m-1)*n), m-1))
C <- diag(m)[, 2:7, drop = FALSE]
known <- list(sigma = 1, weights = 1)
fit <- ridge_glht(Y, X, C, lambda = 1, spectrum = known)
G <- solve(tcrossprod(X))
L <- t(X) %*% G %*% C
P1 <- L %*% solve(crossprod(C, G %*% C)) %*% t(L)
P2 <- diag(n) - t(X) %*% G %*% X
W1 <- Y %*% P1 %*% t(Y) / ncol(C)
W2 <- Y %*% P2 %*% t(Y) / (n-m)
ref <- ridge_root(W1, W2, ncol(C), n-m, lambda = 1, spectrum = known)
near(fit$statistic, ref$statistic)
near(fit$largest.root, max(Re(eigen(W1 %*% solve(W2 + diag(p)), only.values = TRUE)$values)))
# Invariance under invertible reparameterization of the design.
T <- diag(seq_len(m))
near(fit$statistic, ridge_glht(Y, T %*% X, T %*% C, lambda = 1, spectrum = known)$statistic)
for (method in c("bayes_identity", "bayes_covariance", "minimax")) {
    f <- ridge_root(W1, W2, 6, n-m, selection = method,
                    lambda_grid = c(0.3, 1, 3), spectrum = known)
    tab <- f$tuning
    ev <- eigen(W2, symmetric = TRUE, only.values = TRUE)$values
    mm <- vapply(tab$lambda, function(l) mean(1/(ev+l)), numeric(1))
    aa <- 1 - tab$lambda * mm
    rr <- aa/(1 - p/(n-m)*aa)
    expected <- switch(method, bayes_identity = mm, bayes_covariance = rr,
                       minimax = pmin(mm, rr/mean(ev))) / tab$scale
    near(tab$score, expected)
    near(f$parameter, tab$lambda[which.max(expected)])
    stopifnot(f$p.value >= 0, f$p.value <= 1)
}
# Spectral estimation and full public pipelines, including a singular W2.
for (n2 in c(6, 24)) {
    Z <- matrix(rnorm(p*n2), p)
    S <- tcrossprod(Z) / n2
    ev <- pmax(eigen(S, symmetric = TRUE, only.values = TRUE)$values, 0)
    spec <- estimate_spectrum(ev, p/n2, size = 40, grid_size = 60)
    stopifnot(spec$status == 0, spec$loss >= 0, all(spec$weights > 0), all(Im(spec$z_grid) > 0))
    near(sum(spec$weights), 1)
    spec2 <- estimate_spectrum(ev * 5, p/n2, size = 40, grid_size = 60)
    near(spec$sigma * 5, spec2$sigma, 1e-5)
    near(spec$weights, spec2$weights, 1e-5)
    f <- ridge_root(W1, S, 6, n2, selection = "minimax", control = list(size = 40, grid_size = 60))
    stopifnot(is.finite(f$statistic), is.finite(f$p.value))
    # Exercise custom grids and first-order-only estimation.
    custom <- estimate_spectrum(ev, p/n2, order = 1,
        sigma_grid = seq(0, max(ev), length.out = 40), z_grid = spec$z_grid)
    stopifnot(custom$status == 0)
}
errors(ridge_calibration(0, 1, 1, 1))
errors(ridge_calibration(1, 1, 1, 0))
errors(ridge_calibration(1, 1, 1, 1, -1))
errors(estimate_spectrum(c(1, -1), 1))
errors(estimate_spectrum(c(1, 2), 1, z_grid = 1:3))
errors(ridge_root(W1, -W2, 6, 51))
errors(ridge_root(W1, W2, 6, 51, selection = "minimax", lambda = 1))
errors(ridge_glht(Y, rbind(X, X[1, ]), rbind(C, C[1, ])))
errors(ridge_glht(Y, X, cbind(C, C[, 1])))
