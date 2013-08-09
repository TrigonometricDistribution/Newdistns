\name{mog}
\alias{dmog}
\alias{pmog}
\alias{qmog}
\alias{rmog}
\title{Marshall Olkin G distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the Marshall Olkin distribution due to Marshall and Olkin  (1997) specified by the pdf
\deqn{\frac {\displaystyle \beta g (x)}{\displaystyle \left[ 1 - (1 - \beta) G (x) \right]^2}}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf and \eqn{\beta > 0}, the scale parameter.}
\usage{
dmog(x, spec, beta = 1, log = FALSE, ...)
pmog(x, spec, beta = 1, log.p = FALSE, lower.tail = TRUE, ...)
qmog(p, spec, beta = 1, log.p = FALSE, lower.tail = TRUE, ...)
rmog(n, spec, beta = 1, ...)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{beta}{the value of the scale parameter, must be positive, the default is 1}
  \item{spec}{a character string specifying the distribution of G and g (for example, "norm" if G and g correspond to the standard normal).}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
  \item{...}{other parameters}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random numbers generated.}
\references{S. Nadarajah, Newdistns: An R Package for new families of distributions, submitted

A. W. Marshall, I. Olkin, A new method for adding a parameter to a family of distributions with application to the exponential and Weibull families, Biometrika 84 (1997) 641-652}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dmog(x,"exp",beta=1)
pmog(x,"exp",beta=1)
qmog(x,"exp",beta=1)
rmog(10,"exp",beta=1)}
\keyword{Marshall Olkin G distribution}
