\name{gammag1}
\alias{dgammag1}
\alias{pgammag1}
\alias{qgammag1}
\alias{rgammag1}
\title{Gamma G distribution due to Zografos and Balakrishnan (2009)}
\description{Computes the pdf, cdf, quantile and random numbers of the gamma G distribution due to Zografos and Balakrishnan (2009) specified by the pdf
\deqn{f (x) = \frac{\displaystyle 1}{\displaystyle \Gamma(a)} g (x) \left\{ -\log \left[ 1 - G(x) \right] \right\}^{a-1}}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, and \eqn{a > 0}, the shape parameter.}
\usage{
dgammag1(x, spec, a = 1, log = FALSE, ...)
pgammag1(x, spec, a = 1, log.p = FALSE, lower.tail = TRUE, ...)
qgammag1(p, spec, a = 1, log.p = FALSE, lower.tail = TRUE, ...)
rgammag1(n, spec, a = 1, ...)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{a}{the value of the shape parameter, must be positive, the default is 1}
  \item{spec}{a character string specifying the distribution of G and g (for example, "norm" if G and g correspond to the standard normal).}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
  \item{...}{other parameters}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random numbers generated.}
\references{S. Nadarajah, Newdistns: An R Package for new families of distributions, submitted

K. Zografos, N. Balakrishnan, On families of beta- and generalized gamma-generated distributions and associated inference, Statistical Methodology 6 (2009) 344-362}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dgammag1(x,"exp",a=1)
pgammag1(x,"exp",a=1)
qgammag1(x,"exp",a=1)
rgammag1(10,"exp",a=1)}
\keyword{Gamma G distribution due to Zografos and Balakrishnan (2009)}
