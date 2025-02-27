#' One Sample Covariance Test by Srivastava, Yanagihara, and Kubokawa (2014)
#'
#' Given data, it performs 1-sample test for Covariance where
#' the null hypothesis is
#' \deqn{H_0 : \Sigma_n = \Sigma_0}
#' where \eqn{\Sigma_n} is the covariance of data model and \eqn{\Sigma_0} is a
#' hypothesized covariance based on a procedure proposed by Srivastava, Yanagihara, and Kubokawa (2014).
#'
#' @param data an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param Sigma0 a \eqn{(p\times p)} given covariance matrix.
#' @param alpha level of significance.
#' @return a named list containing \describe{
#' \item{statistic}{a test statistic value.}
#' \item{threshold}{rejection criterion to be compared against test statistic.}
#' \item{reject}{a logical; \code{TRUE} to reject null hypothesis, \code{FALSE} otherwise.}
#' }
#' @examples
#' \dontrun{
#' ## generate data from multivariate normal with trivial covariance.
#' p = 5;n=10
#' data = matrix(rnorm(n*p), ncol=p)
#'  alpha=0.05
#' Sigma0=diag(ncol(data))
#'  run the test
#' syk(data, Sigma0, alpha)
#' }
syk <- function(data, Sigma0, alpha){
  n = nrow(data)  # Number of observations (rows in the data matrix)
  p = ncol(data)  # Number of variables (columns in the data matrix)
  # Adjust the matrix
  scaler = 1 / sqrt(Sigma0)  # Scaling factors based on the diagonal of Sigma0
  X.centered = scale(data, center = TRUE, scale = FALSE)  # Center the data
  X.adjusted = sweep(X.centered, 2, scaler, "*")  # Scale the centered data
  # Main computation
  z.val = qnorm(1 - alpha)  # Z-value from the standard normal distribution
  bar.X = matrix(colMeans(data), nrow = p, ncol = 1)  # Column means of the data
  Sn = cov(data)  # Sample covariance matrix
  Yn = data - matrix(c(bar.X), nrow = n, ncol = p, byrow = TRUE)  # Centered data
  a1.hat = sum(diag(Sn)) / p  # First term of the test statistic
  a2.hat = ((n-1)*(n-2)*sum(diag((n-1)^2*Sn%*%Sn)) - n*(n-1)*sum((rowSums(Yn^2))^2) + (sum(diag((n-1)*Sn)))^2) / (p*n*(n-1)*(n-2)*(n-3))  # Second term of the test statistic
  T2 = (n-1)/2 * (a2.hat - 2*a1.hat + 1)  # Test statistic
  # Post-adjusting
  statistic = T2  # Test statistic
  threshold = z.val  # Threshold
  reject = (statistic > threshold)  # Determine whether to reject the null hypothesis
  # Return results
  return(list(statistic=statistic, threshold=threshold, reject=reject))
}
