#'Two-Sample Covariance Test by Yu, Li and Xue (2022)
#'Given two sets of data matrices X and Y, where X is an n1 rows and p cols matrix and Y is an n2 rows and p cols matrix, we conduct hypothesis testing of the covariance matrix between two samples. The null hypothesis is:
#'\deqn{H_0 : \Sigma_1 = \Sigma_2}
#'\eqn{\Sigma_1} and \eqn{\Sigma_2} are the sample covariance matrices of X and Y respectively. This test method is based on the test method proposed by Yu, Li and Xue (2022). When the pval value is less than the significance coefficient (generally 0.05), the null hypothesis is rejected.
#' @param X A matrix of n1 by p.
#' @param Y A matrix of n2 by p.
#'
#' @return
#'\item{stat}{a test statistic value.}
#'\item{pval}{a test p_value.}
#'
#' @examples
#' \dontrun{
#' ## generate data from multivariate normal with trivial covariance.
#'##generate X and Y.
#'p= 500;  n1 = 100; n2 = 150
#'X=matrix(rnorm(n1*p), ncol=p)
#'Y=matrix(rnorm(n2*p), ncol=p)
#'## run test
#'PEF(X,Y)
#' }
PEF <- function (X, Y)
{
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
    warning(paste0("These methods are designed for high-dimensional data!\n                   The data dimension (p=",
                   p, ") is small in the input data,\n                   which may results in an inflated Type-I error rate."))
  }
  output_cov_LC = LC(X, Y)
  output_cov_CLX = CLX(X, Y)
  pval_cov_LC = output_cov_LC$pval
  pval_cov_CLX = output_cov_CLX$pval
  stat_cov_Fisher = -2 * log(pval_cov_LC) - 2 * log(pval_cov_CLX)
  pval_cov_Fisher = 1 - pchisq(stat_cov_Fisher, df = 4)
  return(list(stat = stat_cov_Fisher, pval = pval_cov_Fisher))
}
