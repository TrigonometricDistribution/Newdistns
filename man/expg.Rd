\name{expg}
\alias{dexpg}
\alias{pexpg}
\alias{qexpg}
\alias{rexpg}
\title{Exponentiated G distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the exponentiated G distribution due to Gupta et al. (1998) specified by the pdf
\deqn{f (x) = a g (x) G^{a - 1} (x)}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf and \eqn{a > 0}, the shape parameter.}
\usage{
dexpg(x, spec, a = 1, log = FALSE, ...)
pexpg(x, spec, a = 1, log.p = FALSE, lower.tail = TRUE, ...)
qexpg(p, spec, a = 1, log.p = FALSE, lower.tail = TRUE, ...)
rexpg(n, spec, a = 1, ...)
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

R. C. Gupta, P. L. Gupta, R. D. Gupta, Modeling failure time data by Lehman alternatives, Communications in Statistics---Theory and Methods 27 (1998) 887-904}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dexpg(x,"exp",a=1)
pexpg(x,"exp",a=1)
qexpg(x,"exp",a=1)
rexpg(10,"exp",a=1)}
\keyword{Exponentiated G distribution}
