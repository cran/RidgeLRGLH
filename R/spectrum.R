# Scalar validation shared by the public entry points.
.positive <- function(x, name) {
    if (!is.numeric(x) || is.complex(x) || length(x) != 1L ||
        !is.finite(x) || x <= 0) stop(name, " must be a finite positive number.", call. = FALSE)
    invisible(x)
}
.spectrum <- function(sigma, weights) {
    if (!is.numeric(sigma) || is.complex(sigma) || !length(sigma) ||
        any(!is.finite(sigma)) || any(sigma < 0) || max(sigma) <= 0)
        stop("sigma must be finite, nonnegative, and not all zero.", call. = FALSE)
    if (!is.numeric(weights) || is.complex(weights) || length(weights) != length(sigma) ||
        any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0)
        stop("weights must be finite nonnegative masses matching sigma.", call. = FALSE)
    keep <- weights > 0
    list(sigma = sigma[keep], weights = weights[keep] / sum(weights[keep]))
}
.companion <- function(z, ev, gamma, derivative = FALSE) {
    if (derivative) gamma * mean(1 / (ev - z)^2) - (gamma - 1) / z^2
    else gamma * mean(1 / (ev - z)) + (gamma - 1) / z
}
# Invert the companion transform on each exterior real branch, then continue
# into the upper half-plane to the prescribed transform grid (Supplement S.2).
.transform_grid <- function(ev, gamma, lambda, size) {
    left <- .companion(1.05, ev, gamma)
    right <- .companion(-lambda, ev, gamma)
    targets <- seq(left, right, length.out = size) + 0.01i
    vapply(targets, function(target) {
        re <- Re(target)
        if (re > 0) {
            lo <- -max(1, 2 / re)
            while (.companion(lo, ev, gamma) > re) lo <- lo * 2
            # The positive branch exists down to minus infinity. Our targets
            # do not exceed phi(-lambda).
            z <- stats::uniroot(function(x) .companion(x, ev, gamma) - re,
                                c(lo, -lambda), tol = 1e-11)$root
        } else if (re < 0) {
            hi <- max(2, -2 / re)
            while (.companion(hi, ev, gamma) < re) hi <- hi * 2
            z <- stats::uniroot(function(x) .companion(x, ev, gamma) - re,
                                c(1.05, hi), tol = 1e-11)$root
        } else z <- 1i / Im(target)
        z <- as.complex(z)
        for (k in seq_len(100L)) {
            residual <- .companion(z, ev, gamma) - target
            if (Mod(residual) < 1e-10) break
            step <- residual / .companion(z, ev, gamma, TRUE)
            rate <- 1
            repeat {
                candidate <- z - rate * step
                if (is.finite(Mod(candidate)) && Im(candidate) > 0 &&
                    Mod(.companion(candidate, ev, gamma) - target) < Mod(residual)) break
                rate <- rate / 2
                if (rate < 2^-40) stop("Transform grid inversion failed; supply z_grid.", call. = FALSE)
            }
            z <- candidate
        }
        if (Mod(.companion(z, ev, gamma) - target) > 1e-8)
            stop("Transform grid inversion did not converge; supply z_grid.", call. = FALSE)
        z
    }, complex(1))
}

estimate_spectrum <- function(eigenvalues, gamma2, lambda = NULL,
                              size = 500L, grid_size = 500L,
                              sigma_grid = NULL, z_grid = NULL,
                              order = 2L, truncation = 2) {
    if (!is.numeric(eigenvalues) || is.complex(eigenvalues) || !length(eigenvalues) ||
        any(!is.finite(eigenvalues)) || any(eigenvalues < 0) || max(eigenvalues) <= 0)
        stop("eigenvalues must be finite, nonnegative, and not all zero.", call. = FALSE)
    .positive(gamma2, "gamma2")
    if (is.null(lambda)) lambda <- 0.5 * mean(eigenvalues)
    .positive(lambda, "lambda")
    for (nm in c("size", "grid_size")) {
        val <- get(nm)
        .positive(val, nm)
        if (val < 3 || val != floor(val)) stop(nm, " must be an integer >= 3.", call. = FALSE)
    }
    if (length(order) != 1L || is.na(order) || !order %in% 1:2)
        stop("order must be 1 or 2.", call. = FALSE)
    .positive(truncation, "truncation")
    scale <- max(eigenvalues)
    ev <- eigenvalues / scale
    lam <- lambda / scale
    if (is.null(sigma_grid)) {
        upper <- min(1, -1 / .companion(1.025, ev, gamma2))
        sigma <- seq(max(0, min(ev)), upper, length.out = grid_size)
    } else {
        .spectrum(sigma_grid, rep(1, length(sigma_grid)))
        sigma <- sort(unique(sigma_grid / scale))
    }
    if (is.null(z_grid)) z <- .transform_grid(ev, gamma2, lam, size)
    else {
        if (!is.numeric(z_grid) && !is.complex(z_grid)) stop("z_grid must be numeric or complex.")
        if (!length(z_grid) || any(!is.finite(z_grid)) || any(Im(z_grid) <= 0))
            stop("z_grid must lie in the open upper half-plane.", call. = FALSE)
        z <- z_grid / scale
    }
    phi <- vapply(z, .companion, complex(1), ev = ev, gamma = gamma2)
    phip <- vapply(z, .companion, complex(1), ev = ev, gamma = gamma2, derivative = TRUE)
    Q <- list((z + 1 / phi) / (lam * gamma2))
    if (order == 2L) Q[[2]] <- (1 / phi^2 - 1 / phip) / (gamma2 * lam^2)
    H <- outer(lam * phi, sigma, function(a, b) b / (a * b + lam))
    E <- do.call(rbind, lapply(seq_len(order), function(j) H^j / Mod(Q[[j]])))
    target <- unlist(lapply(Q, function(q) q / Mod(q)))
    A <- rbind(Re(E), Im(E))
    b <- c(Re(target), Im(target))
    if (any(!is.finite(A)) || any(!is.finite(b))) stop("Non-finite spectral constraints.", call. = FALSE)
    K <- length(sigma)
    fit <- lpSolve::lp("min", c(1, rep(0, K)),
        rbind(cbind(-1, A), cbind(-1, -A), c(0, rep(1, K))),
        c(rep("<=", 2 * nrow(A)), "="), c(b, -b, 1))
    if (fit$status != 0L) stop("Spectral linear program failed (status ", fit$status, ").", call. = FALSE)
    weights <- fit$solution[-1L]
    keep <- weights > 10^(-truncation) / K
    if (!any(keep)) stop("Truncation removed all spectral masses.", call. = FALSE)
    list(sigma = sigma[keep] * scale, weights = weights[keep] / sum(weights[keep]),
         loss = fit$objval, status = fit$status, lambda = lambda,
         z_grid = z * scale, order = order)
}

ridge_calibration <- function(sigma, gamma1, gamma2, lambda,
                               weights = rep(1 / length(sigma), length(sigma))) {
    for (nm in c("gamma1", "gamma2", "lambda")) .positive(get(nm), nm)
    spec <- .spectrum(sigma, weights)
    # Normalize the population spectrum for scale-equivariant root tolerances.
    unit <- max(spec$sigma)
    sig <- spec$sigma / unit
    lam <- lambda / unit
    w <- spec$weights
    H <- function(h, j) sum(w * (sig / (lam - sig * h))^j)
    x <- function(h) h + 1 / (1 + gamma2 * H(h, 1))
    xp <- function(h) 1 - gamma2 * H(h, 2) / (1 + gamma2 * H(h, 1))^2
    pole <- lam / max(sig)
    hi <- pole - max(pole * 1e-10, .Machine$double.eps * pole)
    lo <- -1
    while (x(lo) > 0) lo <- lo * 2
    start <- stats::uniroot(x, c(lo, hi), tol = 1e-12)$root
    top_weight <- sum(w[sig == max(sig)])
    if (gamma2 * top_weight < 1) {
        edge <- stats::uniroot(xp, c(start, hi), tol = 1e-12)$root
        rho <- x(edge)
        end <- edge
    } else {
        rho <- pole
        end <- hi
    }
    objective <- function(h) gamma1 * x(h)^2 * H(h, 2) - xp(h)
    h <- stats::uniroot(objective, c(start, end), tol = 1e-12)$root
    beta <- x(h)
    s <- H(h, 1)
    H2 <- H(h, 2)
    H3 <- H(h, 3)
    d <- xp(h)
    xpp <- -2 * gamma2 * H3 / (1 + gamma2 * s)^2 +
        2 * gamma2^2 * H2^2 / (1 + gamma2 * s)^3
    sp <- H2 / d
    spp <- 2 * H3 / d^2 - H2 * xpp / d^3
    out <- c(center = 1 / beta + gamma1 * s,
             scale = (gamma1^3 * spp / 2 + gamma1^2 / beta^3)^(1 / 3),
             beta = beta, rho = rho, s = s, s_prime = sp, s_double_prime = spp)
    if (any(!is.finite(out)) || out["scale"] <= 0 || d <= 0)
        stop("Calibration failed on the regular spectral branch.", call. = FALSE)
    out
}
