\name{weibullg}
\alias{dweibullg}
\alias{pweibullg}
\alias{qweibullg}
\alias{rweibullg}
\title{Weibull G distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the Weibull G distribution due to Alzaatreh et al. (2013) specified by the pdf
\deqn{f (x) = \frac {c}{\beta} \frac {g (x)}{1 - G (x)} \left\{ -\frac {\log \left[ 1 - G (x) \right]}{\beta} \right\}^{c - 1} \exp \left\{ -\left[ -\frac {\log \left[ 1 - G (x) \right]}{\beta} \right]^c \right\}}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, \eqn{\beta > 0}, the scale parameter and \eqn{c > 0}, the shape parameter.}
\usage{
dweibullg(x, spec, beta = 1, c = 1, log = FALSE, ...)
pweibullg(x, spec, beta = 1, c = 1, log.p = FALSE, lower.tail = TRUE, ...)
qweibullg(p, spec, beta = 1, c = 1, log.p = FALSE, lower.tail = TRUE, ...)
rweibullg(n, spec, beta = 1, c = 1, ...)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{beta}{the value of the scale parameter, must be positive, the default is 1}
  \item{c}{the value of the shape parameter, must be positive, the default is 1}
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
dweibullg(x,"exp",beta=1,c=1)
pweibullg(x,"exp",beta=1,c=1)
qweibullg(x,"exp",beta=1,c=1)
rweibullg(10,"exp",beta=1,c=1)}
\keyword{Weibull G distribution}
