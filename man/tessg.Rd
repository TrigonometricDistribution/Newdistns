\name{tessg}
\alias{dtessg}
\alias{ptessg}
\alias{qtessg}
\alias{rtessg}
\title{Truncated-exponential skew-symmetric G distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the truncated-exponential skew-symmetric G distribution due to Nadarajah et al. (2013) specified by the pdf
\deqn{f (x) = \frac {\lambda}{1 - \exp (-\lambda)} g(x) \exp \left\{ -\lambda G(x) \right\}}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, and \eqn{-\infty < \lambda < \infty}, the skewness parameter.}
\usage{
dtessg(x, spec, lambda = 1, log = FALSE, ...)
ptessg(x, spec, lambda = 1, log.p = FALSE, lower.tail = TRUE, ...)
qtessg(p, spec, lambda = 1, log.p = FALSE, lower.tail = TRUE, ...)
rtessg(n, spec, lambda = 1, ...)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{lambda}{the value of skewness parameter, can be any real value, the default is 1}
  \item{spec}{a character string specifying the distribution of G and g (for example, "norm" if G and g correspond to the standard normal).}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
  \item{...}{other parameters}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random numbers generated.}
\references{S. Nadarajah, Newdistns: An R Package for new families of distributions, submitted

S. Nadarajah, V. Nassiri, A. Mohammadpour, Truncated-exponential skew-symmetric distributions, Statistics, to appear}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dtessg(x,"exp",lambda=1)
ptessg(x,"exp",lambda=1)
qtessg(x,"exp",lambda=1)
rtessg(10,"exp",lambda=1)}
\keyword{Truncated-exponential skew-symmetric G distribution}
