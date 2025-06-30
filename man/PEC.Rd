\name{PEC}
\alias{PEC}
\title{Two-Sample Covariance Test by Yu, Li and Xue (2022)}
\usage{
PEC(X,Y)
}
\arguments{
\item{X}{A matrix of n1 by p}
\item{Y}{A matrix of n2 by p}
}
\value{
\item{stat}{a test statistic value.}
\item{pval}{a test p_value.}
}
\references{Yu, X., Li, D., and Xue, L. (2022). Fisher's combined probability test for high-dimensional covariance matrices. Journal of the American Statistical Association, (in press):1-14.}
\description{
Given two sets of data matrices X and Y, where X is an n1 rows and p cols matrix and Y is an n2 rows and p cols matrix,, we conduct hypothesis testing of the covariance matrix between two samples. The null hypothesis is:
\deqn{H_0 : \Sigma_1 = \Sigma_2}
\eqn{\Sigma_1} and \eqn{\Sigma_2} are the sample covariance matrices of X and Y respectively. This test method is based on the test method proposed by Yu, Li and Xue (2022). When the pval value is less than the significance coefficient (generally 0.05), the null hypothesis is rejected.
}
\examples{
## generate X and Y.
p= 500;  n1 = 100; n2 = 150
X=matrix(rnorm(n1*p), ncol=p)
Y=matrix(rnorm(n2*p), ncol=p)
## run test
PEC(X,Y)
}

