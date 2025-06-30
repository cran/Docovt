#'Two-Sample Covariance Test by Li and Chen (2012)
#'Given two sets of data matrices X and Y, where X is an n1 rows and p cols matrix and Y is an n2 rows and p cols matrix, we conduct hypothesis testing of the covariance matrix between two samples. The null hypothesis is:
#'\deqn{H_0 : \Sigma_1 = \Sigma_2}
#'\eqn{\Sigma_1} and \eqn{\Sigma_2} are the sample covariance matrices of X and Y respectively. This test method is based on the test method proposed by Li and Chen (2012). When the pval value is less than the significance coefficient (generally 0.05), the null hypothesis is rejected.
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
#'LC(X,Y)
#' }
LC <- function(X, Y) {
  tr <- function(A) sum(diag(A))
  A_2 <- function(X) {
    n <- nrow(X)
    p <- ncol(X)
    term1 <- (n-1) * sum(colMeans(X)^2)
    term2 <- 2 * sum(colMeans(X) * colSums(X))
    term3 <- sum(X^2)
    (term1 - term2 + term3) / (n*(n-1))
  }

  A_3 <- function(X) {
    n <- nrow(X)
    p <- ncol(X)
    (sum(X)^2 - sum(X^2)) / (n*(n-1)*(n-2))
  }

  C_23 <- function(X, Y) {
    n1 <- nrow(X)
    n2 <- nrow(Y)
    -sum(X) * sum(Y) / (n1*n2*(n1+n2-1))
  }

  C_4 <- function(X, Y) {
    n1 <- nrow(X)
    n2 <- nrow(Y)
    (sum(X)*sum(Y)) / (n1*n2*(n1+n2-1)*(n1+n2-2))
  }
  X = X
  Y = Y
  n1 = nrow(X)
  n2 = nrow(Y)
  p = ncol(X)
  p2 = ncol(Y)

  if (p != p2) {
    stop("The data dimensions ncol(dataX) and ncol(dataY) do not match!")
  }

  if (p <= 30) {
    warning(paste0("These methods are designed for high-dimensional data!\n",
                   "The data dimension (p=", p, ") is small in the input data,\n",
                   "which may result in an inflated Type-I error rate."))
  }
  XXT = X %*% t(X)
  YYT = Y %*% t(Y)
  XYT = X %*% t(Y)

  A.1 = (sum((XXT)^2) - tr((XXT)^2))/(n1*(n1-1))
  B.1 = (sum((YYT)^2) - tr((YYT)^2))/(n2*(n2-1))

  A.2_val = A_2(X)
  A.3_val = A_3(X)
  B.2_val = A_2(Y)
  B.3_val = A_3(Y)

  An1 = A.1 - A.2_val + A.3_val
  Bn2 = B.1 - B.2_val + B.3_val

  C.1 = sum((XYT)^2)/(n1*n2)
  C.2_val = C_23(X, Y)
  C.3_val = C_23(Y, X)
  C.4_val = C_4(X, Y)

  Cn = C.1 + C.2_val + C.3_val + C.4_val
  T.value = An1 + Bn2 - 2*Cn
  sd.est = (2*An1)/n2 + (2*Bn2)/n1
  stat_cov_LC = T.value/sd.est
  pval_cov_LC = 1 - pnorm(stat_cov_LC)

  return(list(stat = stat_cov_LC, pval = pval_cov_LC))
}
