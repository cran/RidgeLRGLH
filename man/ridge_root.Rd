\name{ridge_root}
\alias{ridge_root}
\title{Ridge-Regularized Largest Root Test for Two Covariance Matrices}
\description{Tests using the largest eigenvalue of W1 times the inverse of W2 plus a positive ridge, with asymptotic real Tracy-Widom calibration.}
\usage{
ridge_root(W1, W2, n1, n2, lambda = NULL,
  selection = c("fixed", "bayes_identity", "bayes_covariance", "minimax"),
  lambda_grid = NULL, spectrum = NULL, control = list())
}
\arguments{
\item{W1, W2}{Real symmetric positive semidefinite matrices of equal dimensions, normalized by n1 and n2 respectively. W2 must have positive trace. Under the null they must satisfy the common-covariance double Wishart structure and assumptions of the article.}
\item{n1, n2}{Positive numerator and denominator degrees of freedom.}
\item{lambda}{Positive fixed ridge, defaulting to the mean eigenvalue of W2. Supply only for fixed selection.}
\item{selection}{Fixed ridge, Bayes identity prior, Bayes covariance prior, or minimax over nonnegative mixtures of these priors with equal total signal strength.}
\item{lambda_grid}{Positive candidate ridges for adaptive selection. Defaults to ten equally spaced values from \eqn{\min(0.1p/n_2,1)\operatorname{tr}(W_2)/p} to \eqn{5\operatorname{tr}(W_2)/p}. Ties select the smallest candidate.}
\item{spectrum}{Optional list with population mass locations \code{sigma} and \code{weights}. If omitted, estimates the spectrum from W2. Known spectra may be supplied for oracle experiments.}
\item{control}{Named list of options passed to \code{estimate_spectrum}, excluding \code{eigenvalues} and \code{gamma2}. The \code{lambda} option here controls spectral fitting only.}
}
\details{Bayes criteria are the identity or covariance signal estimates divided by the scaling constant. The minimax criterion is the smaller of the identity signal estimate and the trace-normalized covariance signal estimate, divided by that same scaling constant. This solves the minimization over nonnegative mixtures exactly at an endpoint. All adaptive choices use a finite candidate set. Their asymptotic validity requires the separation and consistency conditions in Section 4 of the article; choosing a ridge from the observed data does not provide a finite-sample error guarantee. The returned p-value is for the selected test and does not adjust for reporting multiple tests or searching multiple hypotheses. Inputs must be real; complex-valued Wishart models are not supported.}
\value{An object of class \code{htest}, with statistic, selected ridge in \code{parameter}, asymptotic upper-tail \code{p.value}, and standard test metadata. Additional fields include \code{largest.root}, \code{center}, \code{scale}, \code{n1}, \code{n2}, \code{dimension}, \code{selection}, \code{tuning} (candidate ridges, constants and scores), \code{spectrum}, and selected \code{calibration}. Fixed selection has an undefined (\code{NA}) tuning score.}
\examples{
set.seed(17)
Z1 <- matrix(rnorm(20 * 30), 20)
Z2 <- matrix(rnorm(20 * 60), 20)
ridge_root(tcrossprod(Z1) / 30, tcrossprod(Z2) / 60,
           n1 = 30, n2 = 60, control = list(size = 40, grid_size = 60))
}
\seealso{\code{\link{ridge_glht}}, \code{\link{estimate_spectrum}}}
\references{Li, H. (2025). Ridge-Regularized Largest Root Test For High-Dimensional General Linear Hypotheses. \doi{10.48550/arXiv.2504.15510}. Accepted by the Annals of Statistics.}
