#'Two-Sample Covariance Test by Cai, Liu and Xia (2013)
#'Given two sets of data matrices X and Y, where X is an n1 rows and p cols matrix and Y is an n2 rows and p cols matrix, we conduct hypothesis testing of the covariance matrix between two samples. The null hypothesis is:
#'\deqn{H_0 : \Sigma_1 = \Sigma_2}
#'\eqn{\Sigma_1} and \eqn{\Sigma_2} are the sample covariance matrices of X and Y respectively. This test method is based on the test method proposed by Cai, Liu and Xia (2013). When the pval value is less than the significance coefficient (generally 0.05), the null hypothesis is rejected.
#'
#' @param X A matrix of n1 by p.
#' @param Y A matrix of n2 by p.
#'
#' @return
#'\item{stat}{a test statistic value.}
#'\item{pval}{a test p_value.}
#'
#'
#' @examples
#' \dontrun{
#' ## generate data from multivariate normal with trivial covariance.
#'##generate X and Y.
#'p= 500;  n1 = 100; n2 = 150
#'X=matrix(rnorm(n1*p), ncol=p)
#'Y=matrix(rnorm(n2*p), ncol=p)
#'## run test
#'CLX(X,Y)
#' }
CLX <- function(X, Y) {
  center <- function(mat) {
    scale(mat, center = TRUE, scale = FALSE)
  }
  F.extreme.cov <- function(x) {
    exp(-1 / sqrt(8 * pi) * exp(-x / 2))
  }

  X = X
  Y = Y
  n1 = nrow(X)
  n2 = nrow(Y)
  p = ncol(X)
  p2 = ncol(Y)
  if (p != p2) {
    stop(" The data dimensions ncol(dataX) and ncol(dataY) do not match!")
  }
  if (p <= 30) {
    warning(paste0("These methods are designed for high-dimensional data!\n",
                   "The data dimension (p=", p, ") is small,\n",
                   "which may result in inflated Type-I error rates."))
  }
  Sigma.hat.1 = ((n1 - 1) / n1) * cov(X)
  Sigma.hat.2 = ((n2 - 1) / n2) * cov(Y)
  X.c = center(X)
  Y.c = center(Y)
  Theta.hat.1 = (1 / n1) * crossprod(X.c^2) -
    (2 / n1) * Sigma.hat.1 * crossprod(X.c) +
    Sigma.hat.1^2
  Theta.hat.2 = (1 / n2) * crossprod(Y.c^2) -
    (2 / n2) * Sigma.hat.2 * crossprod(Y.c) +
    Sigma.hat.2^2
  M.value = max((Sigma.hat.1 - Sigma.hat.2)^2 / (Theta.hat.1 / n1 + Theta.hat.2 / n2))
  stat_cov_CLX = M.value - 4 * log(p) + log(log(p))
  pval_cov_CLX = 1 - F.extreme.cov(stat_cov_CLX)

  return(list(stat = stat_cov_CLX, pval = pval_cov_CLX))
}
