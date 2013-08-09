\name{gbg}
\alias{dgbg}
\alias{pgbg}
\alias{qgbg}
\alias{rgbg}
\title{Generalized beta G distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the generalized beta G distribution due to Alexander et al. (2012) specified by the pdf
\deqn{f (x) = \frac {c}{B (a, b)} g (x) G^{a c - 1} (x) \left[ 1 - G^c (x) \right]^{b - 1}}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, \eqn{a > 0}, the first shape parameter, \eqn{b > 0}, the second shape parameter and \eqn{c > 0}, the third shape parameter.}
\usage{
dgbg(x, spec, a = 1, b = 1, c = 1, log = FALSE, ...)
pgbg(x, spec, a = 1, b = 1, c = 1, log.p = FALSE, lower.tail = TRUE, ...)
qgbg(p, spec, a = 1, b = 1, c = 1, log.p = FALSE, lower.tail = TRUE, ...)
rgbg(n, spec, a = 1, b = 1, c = 1, ...)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{b}{the value of the second shape parameter, must be positive, the default is 1}
  \item{c}{the value of the third shape parameter, must be positive, the default is 1}
  \item{spec}{a character string specifying the distribution of G and g (for example, "norm" if G and g correspond to the standard normal).}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
  \item{...}{other parameters}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random numbers generated.}
\references{S. Nadarajah, Newdistns: An R Package for new families of distributions, submitted

C. Alexander, G. M. Cordeiro, E. M. M. Ortega, Generalized beta-generated distributions, Computational Statistics and Data Analysis 56 (2012) 1880-1897}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dgbg(x,"exp",a=1,b=1,c=1)
pgbg(x,"exp",a=1,b=1,c=1)
qgbg(x,"exp",a=1,b=1,c=1)
rgbg(10,"exp",a=1,b=1,c=1)}
\keyword{Generalized beta G distribution}
