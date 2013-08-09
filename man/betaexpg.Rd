\name{betaexpg}
\alias{dbetaexpg}
\alias{pbetaexpg}
\alias{qbetaexpg}
\alias{rbetaexpg}
\title{Beta exponential G distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the beta exponential G distribution due to Alzaatreh et al. (2013) specified by the pdf
\deqn{f (x) = \frac {\lambda}{B (a, b)} g (x) \left[ 1 - G (x) \right]^{\lambda b - 1}\left\{ 1 - \left[ 1 - G (x) \right]^\lambda \right\}^{a - 1}}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, \eqn{\lambda > 0}, the first shape parameter, \eqn{a >  0}, the second shape parameter, and \eqn{b > 0}, the third shape parameter.}
\usage{
dbetaexpg(x, spec, lambda = 1, a = 1, b = 1, log = FALSE, ...)
pbetaexpg(x, spec, lambda = 1, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE, ...)
qbetaexpg(p, spec, lambda = 1, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE, ...)
rbetaexpg(n, spec, lambda = 1, a = 1, b = 1, ...)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{lambda}{the value of the first parameter, must be positive, the default is 1}
  \item{a}{the value of the second shape parameter, must be positive, the default is 1}
  \item{b}{the value of the third shape parameter, must be positive, the default is 1}
  \item{spec}{a character string specifying the distribution of G and g (for example, "norm" if G and g correspond to the standard normal).}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
  \item{...}{other parameters}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random numbers generated.}
\references{S. Nadarajah, Newdistns: An R Package for new families of distributions, submitted

A. Alzaatreh, C. Lee, F. Famoye, A new method for generating families of continuous distributions, METRON 71 (2013) 63-79}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dbetaexpg(x,"exp",lambda=1,a=1,b=1)
pbetaexpg(x,"exp",lambda=1,a=1,b=1)
qbetaexpg(x,"exp",lambda=1,a=1,b=1)
rbetaexpg(10,"exp",lambda=1,a=1,b=1)}
\keyword{Beta exponential G distribution}
