\name{CLX}
\alias{CLX}
\title{Two-Sample Covariance Test by Cai, Liu and Xia (2013)}
\usage{
CLX(X,Y)
}
\arguments{
\item{X}{A matrix of n1 by p}
\item{Y}{A matrix of n2 by p}
}
\value{
\item{stat}{a test statistic value.}
\item{pval}{a test p_value.}
}
\references{Cai, T. T., Liu, W., and Xia, Y. (2013). Two-sample covariance matrix testing and support recovery in high-dimensional and sparse settings. Journal of the American Statistical Association, 108(501):265-277.}
\description{
Given two sets of data matrices X and Y, where X is an n1 rows and p cols matrix and Y is an n2 rows and p cols matrix, we conduct hypothesis testing of the covariance matrix between two samples. The null hypothesis is:
\deqn{H_0 : \Sigma_1 = \Sigma_2}
\eqn{\Sigma_1} and \eqn{\Sigma_2} are the sample covariance matrices of X and Y respectively. This test method is based on the test method proposed by Cai, Liu and Xia (2013). When the pval value is less than the significance coefficient (generally 0.05), the null hypothesis is rejected.
}
\examples{
## generate X and Y.
p= 500;  n1 = 100; n2 = 150
X=matrix(rnorm(n1*p), ncol=p)
Y=matrix(rnorm(n2*p), ncol=p)
## run test
CLX(X,Y)
}

