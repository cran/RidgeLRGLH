.psd <- function(x, name) {
    if (!is.matrix(x) || !is.numeric(x) || is.complex(x) || nrow(x) != ncol(x) ||
        !nrow(x) || any(!is.finite(x))) stop(name, " must be a finite real square matrix.", call. = FALSE)
    tol <- 1e-8 * max(abs(x), .Machine$double.xmin)
    if (max(abs(x - t(x))) > tol) stop(name, " must be symmetric.", call. = FALSE)
    e <- eigen((x + t(x)) / 2, symmetric = TRUE)
    if (min(e$values) < -tol) stop(name, " must be positive semidefinite.", call. = FALSE)
    e$values <- pmax(e$values, 0)
    e
}

ridge_root <- function(W1, W2, n1, n2, lambda = NULL,
                       selection = c("fixed", "bayes_identity", "bayes_covariance", "minimax"),
                       lambda_grid = NULL, spectrum = NULL,
                       control = list()) {
    selection <- match.arg(selection)
    .positive(n1, "n1")
    .positive(n2, "n2")
    .psd(W1, "W1")
    e <- .psd(W2, "W2")
    if (!identical(dim(W1), dim(W2))) stop("W1 and W2 must have the same dimensions.", call. = FALSE)
    p <- nrow(W2)
    gamma1 <- p / n1
    gamma2 <- p / n2
    tau <- mean(e$values)
    if (tau <= 0) stop("W2 must have positive trace.", call. = FALSE)
    if (!is.list(control) || (length(control) && (is.null(names(control)) ||
        any(!names(control) %in% c("lambda", "size", "grid_size", "sigma_grid", "z_grid", "order", "truncation")))))
        stop("control must be a named list of estimate_spectrum options.", call. = FALSE)
    if (selection == "fixed") {
        if (is.null(lambda)) lambda <- tau
        .positive(lambda, "lambda")
        if (!is.null(lambda_grid)) stop("lambda_grid is only used for adaptive selection.", call. = FALSE)
        grid <- lambda
    } else {
        if (!is.null(lambda)) stop("Use lambda_grid for adaptive selection.", call. = FALSE)
        grid <- lambda_grid
        if (is.null(grid)) grid <- seq(min(0.1 * gamma2, 1), 5, length.out = 10L) * tau
        if (!is.numeric(grid) || is.complex(grid) || !length(grid) || any(!is.finite(grid)) || any(grid <= 0))
            stop("lambda_grid must contain finite positive numbers.", call. = FALSE)
        grid <- sort(unique(grid))
    }
    if (is.null(spectrum)) spectrum <- do.call(estimate_spectrum,
        c(list(eigenvalues = e$values, gamma2 = gamma2), control))
    if (!is.list(spectrum) || is.null(spectrum$sigma) || is.null(spectrum$weights))
        stop("spectrum must contain sigma and weights.", call. = FALSE)
    spec <- .spectrum(spectrum$sigma, spectrum$weights)
    calibration <- t(vapply(grid, function(lam) ridge_calibration(spec$sigma,
        gamma1, gamma2, lam, spec$weights), numeric(7)))
    m <- vapply(grid, function(lam) mean(1 / (e$values + lam)), numeric(1))
    a <- 1 - grid * m
    rho2 <- a / (1 - gamma2 * a)
    score <- switch(selection, fixed = rep(NA_real_, length(grid)),
        bayes_identity = m / calibration[, "scale"],
        bayes_covariance = rho2 / calibration[, "scale"],
        minimax = pmin(m, rho2 / tau) / calibration[, "scale"])
    if (selection != "fixed" && any(!is.finite(score))) stop("Non-finite tuning criterion.", call. = FALSE)
    key <- if (selection == "fixed") 1L else which.max(score)
    chosen <- grid[key]
    # A symmetric similarity transform avoids non-symmetric eigensolvers and
    # preserves the nonnegative real generalized eigenvalues.
    invsqrt <- sweep(e$vectors, 2L, sqrt(e$values + chosen), "/")
    F <- crossprod(invsqrt, W1 %*% invsqrt)
    largest <- max(eigen((F + t(F)) / 2, symmetric = TRUE, only.values = TRUE)$values)
    center <- unname(calibration[key, "center"])
    scaling <- unname(calibration[key, "scale"])
    statistic <- p^(2 / 3) * (largest - center) / scaling
    pvalue <- RMTstat::ptw(statistic, beta = 1, lower.tail = FALSE)
    if (!is.finite(pvalue)) stop("Tracy-Widom probability evaluation failed.", call. = FALSE)
    out <- list(statistic = c(TW = statistic), parameter = c(lambda = chosen),
        p.value = pvalue, alternative = "greater",
        method = paste("Ridge-regularized largest root test (", selection, ")", sep = ""),
        data.name = paste(deparse(substitute(W1)), "and", deparse(substitute(W2))),
        largest.root = largest, center = center, scale = scaling,
        n1 = n1, n2 = n2, dimension = p, selection = selection,
        tuning = data.frame(lambda = grid, center = calibration[, "center"],
                            scale = calibration[, "scale"], score = score),
        spectrum = spectrum, calibration = calibration[key, ])
    class(out) <- "htest"
    out
}

ridge_glht <- function(Y, X, C, lambda = NULL,
                       selection = c("fixed", "bayes_identity", "bayes_covariance", "minimax"),
                       lambda_grid = NULL, spectrum = NULL, control = list()) {
    selection <- match.arg(selection)
    for (nm in c("Y", "X", "C")) {
        val <- get(nm)
        if (!is.matrix(val) || !is.numeric(val) || is.complex(val) ||
            any(!is.finite(val)) || !length(val))
            stop(nm, " must be a nonempty finite real matrix.", call. = FALSE)
    }
    m <- nrow(X)
    n <- ncol(X)
    if (ncol(Y) != n || nrow(C) != m || ncol(C) > m)
        stop("Use Y: p by n, X: m by n, C: m by q with q <= m.", call. = FALSE)
    qx <- qr(t(X))
    if (qx$rank != m || m >= n) stop("X must have full row rank and positive residual degrees of freedom.", call. = FALSE)
    if (qr(C)$rank != ncol(C)) stop("C must have full column rank.", call. = FALSE)
    Q <- qr.Q(qx)
    # t(X) = Q R, hence col(t(X) (X X')^-1 C) = col(Q R'^-1 C).
    R <- qr.R(qx)
    Cperm <- C[qx$pivot, , drop = FALSE]
    A <- Q %*% forwardsolve(t(R), Cperm)
    U <- qr.Q(qr(A))
    residual <- Y - (Y %*% Q) %*% t(Q)
    tested <- Y %*% U
    n1 <- ncol(C)
    n2 <- n - m
    result <- ridge_root(tcrossprod(tested) / n1, tcrossprod(residual) / n2,
        n1, n2, lambda = lambda, selection = selection, lambda_grid = lambda_grid,
        spectrum = spectrum, control = control)
    result$data.name <- paste(deparse(substitute(Y)), "with hypothesis B C = 0")
    result
}
