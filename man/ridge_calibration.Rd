\name{ridge_calibration}
\alias{ridge_calibration}
\title{Compute Tracy-Widom Centering and Scaling Constants}
\description{Computes the spectral edge and calibration constants for a specified population covariance spectrum.}
\usage{
ridge_calibration(sigma, gamma1, gamma2, lambda,
  weights = rep(1 / length(sigma), length(sigma)))
}
\arguments{
\item{sigma}{Nonnegative population covariance eigenvalues or spectral mass locations, not all zero.}
\item{gamma1}{Dimension divided by numerator degrees of freedom.}
\item{gamma2}{Dimension divided by denominator degrees of freedom.}
\item{lambda}{Positive regularization parameter, in covariance units.}
\item{weights}{Nonnegative spectral masses, normalized internally to sum to one.}
}
\details{Uses the scalar parameterization \eqn{x(h)=h+[1+\gamma_2 H_1(h)]^{-1}} from Section 3. It solves for the regular branch and \eqn{\beta^2s'(\beta)=1/\gamma_1} by bracketed root finding. Derivatives follow by implicit differentiation. This is equivalent to the manuscript's ODE characterization, and avoids evaluating a dense ODE trajectory near the edge. Root tolerances are applied after normalization of the population spectrum. The Tracy-Widom approximation assumes the high-dimensional regime and spectral regularity conditions in the article; it is not an exact finite-sample calibration.}
\value{A named numeric vector: \code{center}, \code{scale}, \code{beta}, \code{rho}, \code{s}, \code{s_prime}, and \code{s_double_prime}, with derivatives evaluated at beta. The normalized statistic is \eqn{p^{2/3}(l_{max}-center)/scale}.}
\examples{
ridge_calibration(sigma = 1, gamma1 = 2, gamma2 = 0.5, lambda = 1)
}
\references{Li, H. (2025). Ridge-Regularized Largest Root Test For High-Dimensional General Linear Hypotheses. \doi{10.48550/arXiv.2504.15510}. Accepted by the Annals of Statistics.}
