#' Two Sample Covariance Test by Cai and Ma (2013)
#'
#' Given two sets of data, it performs 2-sample test for equality of covariance matrices where
#' the null hypothesis is
#' \deqn{H_0 : \Sigma_1 = \Sigma_2}
#' where \eqn{\Sigma_1} and \eqn{\Sigma_2} represent true (unknown) covariance
#' for each dataset based on a procedure proposed by Cai and Ma (2013).
#' If \code{statistic} \eqn{>} \code{threshold}, it rejects null hypothesis.
#'
#' @param X an \eqn{(m\times p)}  matrix where each row is an observation from the first dataset.
#' @param Y an \eqn{(n\times p)} matrix where each row is an observation from the second dataset.
#' @param alpha level of significance.
#'
#' @return a named list containing \describe{
#' \item{statistic}{a test statistic value.}
#' \item{threshold}{rejection criterion to be compared against test statistic.}
#' \item{reject}{a logical; \code{TRUE} to reject null hypothesis, \code{FALSE} otherwise.}
#' }
#'
#' @examples
#' ## generate 2 datasets from multivariate normal with identical covariance.
#' p= 5;  n1 = 100; n2 = 150; alpha=0.05
#' X=data1 = matrix(rnorm(n1*p), ncol=p)
#' Y=data2 = matrix(rnorm(n2*p), ncol=p)
#'
#' # run test
#'cmtwo(X, Y, alpha)
#' @export
cmtwo <- function(X, Y, alpha){
  # Parameter setting
  n1 = nrow(X)  # Number of observations in the first sample
  n2 = nrow(Y)  # Number of observations in the second sample
  p  = ncol(X)   # Number of variables
  # Elementary computation: sample covariance with new degrees of freedom
  Sigma1.hat = cov(X) * (n1 - 1) / n1  # Sample covariance of the first sample
  Sigma2.hat = cov(Y) * (n2 - 1) / n2  # Sample covariance of the second sample
  # Mean estimation
  bar.X = colMeans(X)  # Column means of the first sample
  theta1.hat = matrix(0, nrow = p, ncol = p)  # Initialize matrix for the first sample
  for (i in 1:p) {
    for (j in 1:p) {
     theta1.hat[i, j] = sum((((X[, i] - bar.X[i]) * (X[, j] - bar.X[j])) - Sigma1.hat[i, j])^2) / n1
   }
  }
  bar.Y = colMeans(Y)  # Column means of the second sample
  theta2.hat = matrix(0, nrow = p, ncol = p)  # Initialize matrix for the second sample
  for (i in 1:p) {
    for (j in 1:p) {
     theta2.hat[i, j] = sum((((Y[, i] - bar.Y[i]) * (Y[, j] - bar.Y[j])) - Sigma2.hat[i, j])^2) / n2
   }
  }
  # Statistic and results
  M.mat = (Sigma1.hat - Sigma2.hat)^2 / (theta1.hat / n1 + theta2.hat / n2)
  Mn = max(M.mat)  # Test statistic
  q.alpha = -log(8 * pi) - 2 * log(log((1 / alpha)))
  threshold = q.alpha + 4 * log(p) - log(log(p))
  res = (Mn > threshold)  # Decision
  # Return results
  return(list(statistic = Mn, threshold = threshold, reject = res))
}
