\name{eg}
\alias{deg}
\alias{peg}
\alias{qeg}
\alias{reg}
\title{Exponentiated generalized G distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the exponentiated generalized G distribution due to Cordeiro et al. (2013) specified by the pdf
\deqn{f (x) = a b g (x) \left[ 1 - G (x) \right]^{a - 1} \left\{ 1 - \left[ 1 - G (x) \right]^a \right\}^{b - 1}}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, \eqn{a > 0}, the first shape parameter, and \eqn{b > 0}, the second shape parameter.}
\usage{
deg(x, spec, a = 1, b = 1, log = FALSE, ...)
peg(x, spec, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE, ...)
qeg(p, spec, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE, ...)
reg(n, spec, a = 1, b = 1, ...)
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

G. M. Cordeiro, E. M. M. Ortega, D. C. C. da Cunha, The exponentiated generalized class of distributions, Journal of Data Science 11 (2013) 1-27}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
deg(x,"exp",a=1,b=1)
peg(x,"exp",a=1,b=1)
qeg(x,"exp",a=1,b=1)
reg(10,"exp",a=1,b=1)}
\keyword{Exponentiated generalized G distribution}
