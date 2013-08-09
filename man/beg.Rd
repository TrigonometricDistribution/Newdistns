\name{beg}
\alias{dbeg}
\alias{pbeg}
\alias{qbeg}
\alias{rbeg}
\title{Beta extended G distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the beta extended distribution due to Cordeiro et al. (2012) specified by the pdf
\deqn{f (x) = \frac {\displaystyle \alpha g (x)}{\displaystyle B (a, b)}
\left\{ 1 - \exp \left[ -\alpha G (x) \right] \right\}^{a - 1} \exp \left[ -\alpha b G (x) \right]}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, \eqn{\alpha > 0}, the scale parameter, \eqn{a > 0}, the first shape parameter, and \eqn{b > 0}, the second shape parameter.}
\usage{
dbeg(x, spec, alpha = 1, a = 1, b = 1, log = FALSE, ...)
pbeg(x, spec, alpha = 1, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE, ...)
qbeg(p, spec, alpha = 1, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE, ...)
rbeg(n, spec, alpha = 1, a = 1, b = 1, ...)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{alpha}{the value of the scale parameter, must be positive, the default is 1}
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

G. M. Cordeiro, E. M. M. Ortega, G. Silva, The beta extended Weibull family, Journal of Probability and Statistical Science 10 (2012) 15-40}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dbeg(x,"exp",alpha=1,a=1,b=1)
pbeg(x,"exp",alpha=1,a=1,b=1)
qbeg(x,"exp",alpha=1,a=1,b=1)
rbeg(10,"exp",alpha=1,a=1,b=1)}
\keyword{Beta extended G distribution}
