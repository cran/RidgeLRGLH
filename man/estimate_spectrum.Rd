\name{estimate_spectrum}
\alias{estimate_spectrum}
\title{Estimate a Population Covariance Spectrum by Linear Programming}
\description{Fits nonnegative spectral masses using the complex-valued relative discrepancy constraints in Algorithm 1 of Li (2025).}
\usage{
estimate_spectrum(eigenvalues, gamma2, lambda = NULL,
  size = 500L, grid_size = 500L, sigma_grid = NULL,
  z_grid = NULL, order = 2L, truncation = 2)
}
\arguments{
\item{eigenvalues}{All eigenvalues of the normalized residual covariance matrix, including zeros, in any order.}
\item{gamma2}{Dimension divided by residual degrees of freedom.}
\item{lambda}{Positive ridge used for fitting. Defaults to half the mean eigenvalue. This fitting parameter need not equal the eventual test ridge.}
\item{size}{Number of complex transform grid points, at least 3.}
\item{grid_size}{Number of population spectrum grid points, at least 3.}
\item{sigma_grid}{Optional nonnegative population spectrum grid, in original covariance units.}
\item{z_grid}{Optional complex points in the open upper half-plane, in original covariance units. Overrides the default transform grid.}
\item{order}{1 for first-order constraints; 2 (default) also uses derivative constraints.}
\item{truncation}{Positive exponent: discard weights at most \eqn{10^{-truncation}/K}, where K is the spectrum grid size.}
}
\details{The default grid is obtained by inverting the empirical companion Stieltjes transform. Its real parts run from the transform at 1.05 times the largest sample eigenvalue to the transform at minus the fitting ridge. The imaginary transform component is 0.01 divided by the largest eigenvalue. The default population grid follows Supplement S.2. A bounded Newton iteration with step halving performs complex inversion; failed inversion or optimization raises an error. The returned loss is measured before truncation. Finite grids are approximations; inspect sensitivity to grid size, fitting ridge and constraint order.}
\value{A list containing positive \code{weights}, their \code{sigma} locations, pre-truncation \code{loss}, solver \code{status} (zero on success), fitting \code{lambda}, \code{z_grid}, and constraint \code{order}.}
\examples{
set.seed(12)
Z <- matrix(rnorm(20 * 60), 20)
ev <- eigen(tcrossprod(Z) / 60, symmetric = TRUE, only.values = TRUE)$values
spectrum <- estimate_spectrum(ev, 20 / 60, size = 40, grid_size = 60)
sum(spectrum$weights)
}
\references{Li, H. (2025). Ridge-Regularized Largest Root Test For High-Dimensional General Linear Hypotheses. \doi{10.48550/arXiv.2504.15510}. Accepted by the Annals of Statistics.}
