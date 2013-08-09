\name{eepg}
\alias{deepg}
\alias{peepg}
\alias{qeepg}
\alias{reepg}
\title{Exponentiated exponential Poisson G distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the exponentiated exponential Poisson G distribution due to Ristic and Nadarajah (2013) specified by the pdf
\deqn{f (x) = \displaystyle a \lambda \left\{ 1 - \exp (-\lambda)  \right\}^{-1} g (x) G^{a - 1} (x) \exp \left[ -\lambda G^a (x) \right]}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, \eqn{\lambda > 0}, the scale parameter, and \eqn{a > 0}, the  shape parameter.}
\usage{
deepg(x, spec, lambda = 1, a = 1, log = FALSE, ...)
peepg(x, spec, lambda = 1, a = 1, log.p = FALSE, lower.tail = TRUE, ...)
qeepg(p, spec, lambda = 1, a = 1, log.p = FALSE, lower.tail = TRUE, ...)
reepg(n, spec, lambda = 1, a = 1, ...)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{lambda}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the shape parameter, must be positive, the default is 1}
  \item{spec}{a character string specifying the distribution of G and g (for example, "norm" if G and g correspond to the standard normal).}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
  \item{...}{other parameters}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random numbers generated.}
\references{S. Nadarajah, Newdistns: An R Package for new families of distributions, submitted

M. M. Ristic, S. Nadarajah, A new lifetime distribution, Journal of Statistical Computation and Simulation, doi: 10.1080/00949655.2012.697163}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
deepg(x,"exp",lambda=1,a=1)
peepg(x,"exp",lambda=1,a=1)
qeepg(x,"exp",lambda=1,a=1)
reepg(10,"exp",lambda=1,a=1)}
\keyword{Exponentiated exponential Poisson G distribution}
