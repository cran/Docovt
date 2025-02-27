\name{cmtwo}
\alias{cmtwo}
\title{Two-Sample Covariance Test by Cai and Ma (2013)}
\usage{
cmtwo(X, Y, alpha)
}
\arguments{
\item{X}{an \eqn{(m\times p)}  matrix where each row is an observation from the first dataset.}

\item{Y}{an \eqn{(n\times p)} matrix where each row is an observation from the second dataset.}

\item{alpha}{level of significance.}
}
\value{
a named list containing \describe{
\item{statistic}{a test statistic value.}
\item{threshold}{rejection criterion to be compared against test statistic.}
\item{reject}{a logical; \code{TRUE} to reject null hypothesis, \code{FALSE} otherwise.}
}
}
\description{
Given two sets of data, it performs 2-sample test for equality of covariance matrices where
the null hypothesis is
\deqn{H_0 : \Sigma_1 = \Sigma_2}
where \eqn{\Sigma_1} and \eqn{\Sigma_2} represent true (unknown) covariance
for each dataset based on a procedure proposed by Cai and Ma (2013).
If \code{statistic} \eqn{>} \code{threshold}, it rejects null hypothesis.
}
\examples{
## generate 2 datasets from multivariate normal with identical covariance.
p= 5;  n1 = 100; n2 = 150; alpha=0.05
X=data1 = matrix(rnorm(n1*p), ncol=p)
Y=data2 = matrix(rnorm(n2*p), ncol=p)
 
# run test
cmtwo(X, Y, alpha)
}

