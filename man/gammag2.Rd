\name{gammag2}
\alias{dgammag2}
\alias{pgammag2}
\alias{qgammag2}
\alias{rgammag2}
\title{Gamma G distribution due to Ristic and Balakrishnan (2012)}
\description{Computes the pdf, cdf, quantile and random numbers of the gamma G distribution due to Ristic and Balakrishnan (2012) specified by the pdf
\deqn{f (x) = \frac {\displaystyle 1}{\displaystyle \Gamma(a)} g (x) \left\{ -\log G(x) \right\}^{a-1}}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, and \eqn{a > 0}, the shape parameter.}
\usage{
dgammag2(x, spec, a = 1, log = FALSE, ...)
pgammag2(x, spec, a = 1, log.p = FALSE, lower.tail = TRUE, ...)
qgammag2(p, spec, a = 1, log.p = FALSE, lower.tail = TRUE, ...)
rgammag2(n, spec, a = 1, ...)
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

M. M. Ristic, N. Balakrishnan, The gamma exponentiated exponential distribution, Journal of Statistical Computation and Simulation 82 (2012) 1191-1206}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dgammag2(x,"exp",a=1)
pgammag2(x,"exp",a=1)
qgammag2(x,"exp",a=1)
rgammag2(10,"exp",a=1)}
\keyword{Gamma G distribution due to Ristic and Balakrishnan (2012)}
