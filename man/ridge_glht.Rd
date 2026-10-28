\name{ridge_glht}
\alias{ridge_glht}
\title{Test a General Linear Hypothesis in Multivariate Regression}
\description{Tests the hypothesis BC = 0 in the model Y = BX + E using a ridge-regularized largest root statistic.}
\usage{
ridge_glht(Y, X, C, lambda = NULL,
  selection = c("fixed", "bayes_identity", "bayes_covariance", "minimax"),
  lambda_grid = NULL, spectrum = NULL, control = list())
}
\arguments{
\item{Y}{Real p by n response matrix: variables in rows and observations in columns.}
\item{X}{Real m by n fixed design matrix with full row rank, m < n. Include an intercept explicitly if needed.}
\item{C}{Real m by q full column rank contrast matrix specifying BC = 0.}
\item{lambda}{Positive fixed ridge, defaulting to the mean eigenvalue of W2. Supply only for fixed selection.}
\item{selection}{Fixed ridge, Bayes identity prior, Bayes covariance prior, or minimax over nonnegative mixtures of these priors with equal total signal strength.}
\item{lambda_grid}{Positive candidate ridges for adaptive selection. Defaults to ten equally spaced values from \eqn{\min(0.1p/n_2,1)\operatorname{tr}(W_2)/p} to \eqn{5\operatorname{tr}(W_2)/p}. Ties select the smallest candidate.}
\item{spectrum}{Optional list with population mass locations \code{sigma} and \code{weights}. If omitted, estimates the spectrum from W2. Known spectra may be supplied for oracle experiments.}
\item{control}{Named list of options passed to \code{estimate_spectrum}, excluding \code{eigenvalues} and \code{gamma2}. The \code{lambda} option here controls spectral fitting only.}
}
\details{QR factorizations construct the tested and residual subspaces, without forming observation-by-observation projection matrices. Numerator and denominator degrees of freedom are q and n - m. Nonzero right-hand sides and rank-deficient designs or contrasts are not supported. This is an asymptotic test under the moment, covariance, dimensional and regular-edge conditions in Li (2025). Very small tested ranks or dimensions may give poor finite-sample approximation.}
\value{An object of class \code{htest}, with statistic, selected ridge in \code{parameter}, asymptotic upper-tail \code{p.value}, and standard test metadata. Additional fields include \code{largest.root}, \code{center}, \code{scale}, \code{n1}, \code{n2}, \code{dimension}, \code{selection}, \code{tuning} (candidate ridges, constants and scores), \code{spectrum}, and selected \code{calibration}. Fixed selection has an undefined (\code{NA}) tuning score.}
\examples{
set.seed(8)
Y <- matrix(rnorm(20 * 80), 20)
X <- rbind(1, matrix(rnorm(10 * 80), 10))
C <- diag(11)[, -1, drop = FALSE]
ridge_glht(Y, X, C, selection = "minimax",
           control = list(size = 40, grid_size = 60))
}
\seealso{\code{\link{ridge_root}}}
\references{Li, H. (2025). Ridge-Regularized Largest Root Test For High-Dimensional General Linear Hypotheses. \doi{10.48550/arXiv.2504.15510}. Accepted by the Annals of Statistics.}
