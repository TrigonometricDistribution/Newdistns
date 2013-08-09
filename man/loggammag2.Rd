\name{loggammag2}
\alias{dloggammag2}
\alias{ploggammag2}
\alias{qloggammag2}
\alias{rloggammag2}
\title{Log gamma G II distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the log gamma G II distribution due to Amini et al. (2013) specified by the pdf
\deqn{f (x) = \frac {b^a}{\Gamma (a)} g (x) \left\{ -\log G (x) \right\}^{a - 1} G^{b - 1} (x)}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, \eqn{a > 0}, the first shape parameter, \eqn{b > 0}, and the second shape parameter.}
\usage{
dloggammag2(x, spec, a = 1, b = 1, log = FALSE, ...)
ploggammag2(x, spec, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE, ...)
qloggammag2(p, spec, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE, ...)
rloggammag2(n, spec, a = 1, b = 1,  ...)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{b}{the value of the second shape parameter, must be positive, the default is 1}
  \item{spec}{a character string specifying the distribution of G and g (for example, "norm" if G and g correspond to the standard normal).}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
  \item{...}{other parameters}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random numbers generated.}
\references{S. Nadarajah, Newdistns: An R Package for new families of distributions, submitted

M. Amini, S. M. T. K. MirMostafaee, J. Ahmadi, Log-gamma-generated families of distributions, Statistics, 2013, doi: 10.1080/02331888.2012.748775}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dloggammag2(x,"exp",a=1,b=1)
ploggammag2(x,"exp",a=1,b=1)
qloggammag2(x,"exp",a=1,b=1)
rloggammag2(10,"exp",a=1,b=1)}
\keyword{Log gamma G II distribution}
