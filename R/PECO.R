#'Two-Sample Covariance Test by Yu, Li, Xue and Li(2022)
#'Given two sets of data matrices X and Y, where X is an n1 rows and p cols matrix and Y is an n2 rows and p cols matrix, we conduct hypothesis testing of the covariance matrix between two samples. The null hypothesis is:
#'\deqn{H_0 : \Sigma_1 = \Sigma_2}
#'\eqn{\Sigma_1} and \eqn{\Sigma_2} are the sample covariance matrices of X and Y respectively. This test method is based on the test method proposed by Yu, Li, Xue and Li (2022). When the pval value is less than the significance coefficient (generally 0.05), the null hypothesis is rejected.
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
#'PECO(X,Y)
#' }
PECO <- function (X, Y, delta = NULL)
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
  if (is.null(delta)) {
    delta = 4 * log(log(n1 + n2)) * log(p)
  }
  XTX = t(X) %*% X
  X2TX2 = t(X^2) %*% (X^2)
  YTY = t(Y) %*% Y
  Y2TY2 = t(Y^2) %*% (Y^2)
  Aij1 = (XTX^2 - X2TX2)/(n1 * (n1 - 1))
  Bij1 = (YTY^2 - Y2TY2)/(n2 * (n2 - 1))
  Cij1 = (XTX * YTY)/(n1 * n2)
  xbar = apply(X, 2, mean)
  ybar = apply(Y, 2, mean)
  xbarmat = matrix(xbar, p, n1, byrow = F)
  ybarmat = matrix(ybar, p, n2, byrow = F)
  Aij2 = 2/(n1 * (n1 - 1) * (n1 - 2)) * (n1 * n1 * (XTX) *
                                           (xbar %*% t(xbar)) - n1 * (t(X) * xbarmat) %*% (X^2) -
                                           n1 * t(X^2) %*% (X * t(xbarmat)) - XTX * XTX + 2 * X2TX2)
  Bij2 = 2/(n2 * (n2 - 1) * (n2 - 2)) * (n2 * n2 * (YTY) *
                                           (ybar %*% t(ybar)) - n2 * (t(Y) * ybarmat) %*% (Y^2) -
                                           n2 * t(Y^2) %*% (Y * t(ybarmat)) - YTY * YTY + 2 * Y2TY2)
  Cij2 = -((n1 * n1 * xbar %*% t(xbar) - XTX) * (YTY))/(n1 *
                                                          n2 * (n1 - 1))
  Cij3 = -((n2 * n2 * ybar %*% t(ybar) - YTY) * (XTX))/(n1 *
                                                          n2 * (n2 - 1))
  x2bar = apply(X^2, 2, mean)
  y2bar = apply(Y^2, 2, mean)
  Aij3 = (n1^4 * (xbar %*% t(xbar))^2 + -n1 * n1 * XTX * (xbar %*%
                                                            t(xbar)) + -n1 * n1 * n1 * ((xbar * xbar) %*% t(x2bar)) +
            -n1 * n1 * (xbar %*% t(xbar)) * XTX + -n1 * n1 * n1 *
            (x2bar %*% t(xbar * xbar)) + -n1 * n1 * (xbar %*% t(xbar)) *
            XTX + 2 * n1 * t(X^2) %*% (X * t(xbarmat)) + 2 * n1 *
            (t(X) * xbarmat) %*% (X^2) + n1 * n1 * x2bar %*% t(x2bar) +
            XTX * XTX + -n1 * n1 * XTX * (xbar %*% t(xbar)) + 2 *
            n1 * (t(X) * xbarmat) %*% (X^2) + 2 * n1 * (t(X^2) %*%
                                                          (X * t(xbarmat))) + XTX * XTX - 6 * X2TX2)/(n1 * (n1 -
                                                                                                              1) * (n1 - 2) * (n1 - 3))
  Bij3 = (n2^4 * (ybar %*% t(ybar))^2 + -n2 * n2 * YTY * (ybar %*%
                                                            t(ybar)) + -n2 * n2 * n2 * ((ybar * ybar) %*% t(y2bar)) +
            -n2 * n2 * (ybar %*% t(ybar)) * YTY + -n2 * n2 * n2 *
            (y2bar %*% t(ybar * ybar)) + -n2 * n2 * (ybar %*% t(ybar)) *
            YTY + 2 * n2 * t(Y^2) %*% (Y * t(ybarmat)) + 2 * n2 *
            (t(Y) * ybarmat) %*% (Y^2) + n2 * n2 * y2bar %*% t(y2bar) +
            YTY * YTY + -n2 * n2 * YTY * (ybar %*% t(ybar)) + 2 *
            n2 * (t(Y) * ybarmat) %*% (Y^2) + 2 * n2 * (t(Y^2) %*%
                                                          (Y * t(ybarmat))) + YTY * YTY - 6 * Y2TY2)/(n2 * (n2 -
                                                                                                              1) * (n2 - 2) * (n2 - 3))
  Cij4 = (n1 * n1 * xbar %*% t(xbar) - XTX) * (n2 * n2 * ybar %*%
                                                 t(ybar) - YTY)/(n1 * n2 * (n1 - 1) * (n2 - 1))
  Anij = Aij1 - Aij2 + Aij3
  Bnij = Bij1 - Bij2 + Bij3
  Cnij = Cij1 + Cij2 + Cij3 + Cij4
  Tij = Anij + Bnij - 2 * Cnij
  covx = cov(X)
  covy = cov(Y)
  varx = diag(covx)
  vary = diag(covy)
  varxmat = varx %*% t(varx)
  varymat = vary %*% t(vary)
  varmat = 2 * ((1/n1) * covx^2 + (1/n2) * covy^2 + (1/n1) *
                  varxmat + (1/n2) * varymat)^2
  sdmat = sqrt(varmat)
  stat_normalized = Tij/sdmat
  delta = 4 * log(log(n1 + n2)) * log(p)
  idx = sqrt(2) * stat_normalized + 1 > delta
  J0 = sqrt(p) * sum(stat_normalized[idx])
  stat_cov_LC = LC(X, Y)$stat
  stat_cov_PE = stat_cov_LC + J0
  pval_cov_PE = 1 - pnorm(stat_cov_PE)
  return(list(stat = stat_cov_PE, pval = pval_cov_PE))
}
