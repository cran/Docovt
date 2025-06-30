#' One-Sample Covariance Test by Cai and Ma (2013)
#' Given data, it performs 1-sample test for Covariance where
#' the null hypothesis is
#' \deqn{H_0 : \Sigma_n = \Sigma_0}
#' where \eqn{\Sigma_n} is the covariance of data model and \eqn{\Sigma_0} is a
#' hypothesized covariance based on a procedure proposed by Cai and Ma (2013).
#'
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param Sigma0 a \eqn{(p\times p)} given covariance matrix.
#' @param alpha level of significance.
#'
#' @return a named list containing: \describe{
#' \item{statistic}{a test statistic value.}
#' \item{threshold}{rejection criterion to be compared against test statistic.}
#' \item{reject}{a logical; \code{TRUE} to reject null hypothesis, \code{FALSE} otherwise.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate data from multivariate normal with trivial covariance.
#' p = 5;n=10
#' X=data= matrix(rnorm(n*p), ncol=p)
#'  alpha=0.05
#' Sigma0=diag(ncol(X))
#' cm13(X,Sigma0, alpha)
#' }
cm13<- function(X, Sigma0, alpha){
  n = nrow(X)
  p = ncol(X)
 # Adjust the matrix
  scaler = (1 / sqrt(Sigma0))
  X.centered = scale(X, center=TRUE, scale=FALSE)
  X.adjusted = (matrix(X.centered,nrow=n) %*% scaler)
   # Main computation
    # Define the dot product function
dot <- function(a, b) {
  return(sum(a * b))}
  # Initialize hXiXj
  hXiXj <- 0
  # Iterate to compute hXiXj
  for (j in 1:(n-1)) {
    Xj <- X[j, ]
    for (i in 1:(j-1)) {
      Xi <- X[i, ]
      hXiXj <- hXiXj + (dot(Xi, Xj))^2 - (dot(Xi, Xi) + dot(Xj, Xj)) + p
    }
  }
  hXiXj
# Post-adjusting
  statistic = (hXiXj*2)/(n*(n-1)) # test statistic
  threshold = (qnorm(1-alpha))*2*sqrt((p*(p+1))/(n*(n-1)))
  reject = (statistic > threshold)
    # Return results
 return(list(statistic=statistic, threshold=threshold, reject=reject))
}
