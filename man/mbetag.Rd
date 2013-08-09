\name{mbetag}
\alias{dmbetag}
\alias{pmbetag}
\alias{qmbetag}
\alias{rmbetag}
\title{Modified beta G distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the modified beta G distribution due  to Nadarajah et al. (2013) specified by the pdf
\deqn{f (x) = \frac {\displaystyle \beta^a}{\displaystyle B(a,b)} \frac {\displaystyle g(x)\left[ G(x) \right]^{a - 1} \left[ 1 - G(x) \right]^{b - 1}}{\displaystyle \left[ 1 - (1 - \beta ) G(x) \right]^{a + b}}}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, \eqn{\beta > 0}, the scale parameter, \eqn{a > 0}, the first shape parameter and \eqn{b > 0}, the second shape parameter.}
\usage{
dmbetag(x, spec, beta = 1, a = 1, b = 1, log = FALSE, ...)
pmbetag(x, spec, beta = 1, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE, ...)
qmbetag(p, spec, beta = 1, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE, ...)
rmbetag(n, spec, beta = 1, a = 1, b = 1, ...)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{beta}{the value of the scale parameter, must be positive, the default is 1}
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

S. Nadarajah, M. Teimouri, S. H. Shih, Modified beta distributions, Sankhya, 2013, to appear}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dmbetag(x,"exp",beta=1,a=1,b=1)
pmbetag(x,"exp",beta=1,a=1,b=1)
qmbetag(x,"exp",beta=1,a=1,b=1)
rmbetag(10,"exp",beta=1,a=1,b=1)}
\keyword{Modified beta G distribution}
