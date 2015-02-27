\name{gepg}
\alias{dgepg}
\alias{pgepg}
\alias{qgepg}
\alias{rgepg}
\alias{mgepg}
\title{Geometric Exponential Poisson G Distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the geometric exponential Poisson G distribution due to Nadarajah et al. (2013) specified by the pdf
\deqn{f (x) = \displaystyle \frac {\displaystyle \theta (1 - \eta) \left[ 1 - \exp (-\theta) \right] g (x) \exp \left[ -\theta + \theta G (x) \right]}{\displaystyle \left\{ 1 - \exp (-\theta) - \eta + \eta \exp \left[ -\theta + \theta G (x) \right] \right\}^2}}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, \eqn{\theta > 0}, the first scale parameter, and \eqn{0 < eta < 1}, the  second scale parameter.  Also computes the Cramer-von Misses statistic, Anderson Darling statistic, Kolmogorov Smirnov test statistic and p-value, maximum likelihood estimates, Akaike Information Criterion, Consistent Akaikes Information Criterion, Bayesian Information Criterion, Hannan-Quinn information criterion, standard errors of the maximum likelihood estimates, minimum value of the negative log-likelihood function and convergence status when the distribution is fitted to some data}
\usage{
dgepg(x, spec, theta = 1, eta = 0.5, log = FALSE, ...)
pgepg(x, spec, theta = 1, eta = 0.5, log.p = FALSE, lower.tail = TRUE, ...)
qgepg(p, spec, theta = 1, eta = 0.5, log.p = FALSE, lower.tail = TRUE, ...)
rgepg(n, spec, theta = 1, eta = 0.5, ...)
mgepg(g, data, starts, method = "BFGS")
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of probabilities at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{theta}{the value of first scale parameter, must be positive, the default is 1}
  \item{eta}{the value of second scale parameter, must be in the open unit interval, the default is 0.5}
  \item{spec}{a character string specifying the distribution of G and g (for example, "norm" if G and g correspond to the standard normal).}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
  \item{...}{other parameters}
  \item{g}{same as spec but must be one of chisquare ("chisq"), exponential ("exp"), F ("f"), gamma ("gamma"), lognormal ("lognormal"), Weibull ("weibull"), Burr XII ("burrxii"), Chen ("chen"), Frechet ("frechet"), Gompertz ("gompertz"), linear failure rate ("lfr"), log-logistic ("log-logistic"), Lomax ("lomax") and Rayleigh ("rayleigh").  Each of these distributions has one parameter (\code{r}) or two parameters (\code{r}, \code{s}),  for details  including the density function and parameter specifications see Nadarajah and Rocha (2014)}
  \item{data}{a vector of data values for which the distribution is to be fitted}
  \item{starts}{initial values of \code{(theta, eta, r)} if \code{g} has one parameter or initial values of \code{(theta, eta, r, s)} if \code{g} has two parameters}
  \item{method}{the method for optimizing the log likelihood function.  It can be one of \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"} or \code{"SANN"}.  The default is \code{"BFGS"}.  The details of these methods can be found in the manual pages for  \code{optim}}}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random numbers generated or an object giving the values of Cramer-von Misses statistic, Anderson Darling statistic, Kolmogorov Smirnov test statistic and p-value, maximum likelihood estimates, Akaike Information Criterion, Consistent Akaikes Information Criterion, Bayesian Information Criterion, Hannan-Quinn information criterion, standard errors of the maximum likelihood estimates, minimum value of the negative log-likelihood function and convergence status.}
\references{S. Nadarajah and R. Rocha, Newdistns: An R Package for new families of distributions, submitted

S. Nadarajah, V. G. Cancho, E. M. M. Ortega, The geometric exponential Poisson distribution, Stat Methods Appl 22 (2013) 355-380}
\author{Saralees Nadarajah, Ricardo Rocha}
\examples{x=runif(10,min=0,max=1)
dgepg(x,"exp",theta=1,eta=0.5)
pgepg(x,"exp",theta=1,eta=0.5)
qgepg(x,"exp",theta=1,eta=0.5)
rgepg(10,"exp",theta=1,eta=0.5)
mgepg("exp",rexp(100),starts=c(1,0.5,1),method="BFGS")}
\keyword{Geometric exponential Poisson G distribution}
