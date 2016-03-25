\name{eepg}
\alias{deepg}
\alias{peepg}
\alias{qeepg}
\alias{reepg}
\alias{meepg}
\title{Exponentiated Exponential Poisson G Distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the exponentiated exponential Poisson G distribution due to Ristic and Nadarajah (2013) specified by the pdf
\deqn{f (x) = \displaystyle a \lambda \left\{ 1 - \exp (-\lambda)  \right\}^{-1} g (x) G^{a - 1} (x) \exp \left[ -\lambda G^a (x) \right]}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, \eqn{\lambda > 0}, the scale parameter, and \eqn{a > 0}, the  shape parameter.  Also computes the Cramer-von Misses statistic, Anderson Darling statistic, Kolmogorov Smirnov test statistic and p-value, maximum likelihood estimates, Akaike Information Criterion, Consistent Akaikes Information Criterion, Bayesian Information Criterion, Hannan-Quinn information criterion, standard errors of the maximum likelihood estimates, minimum value of the negative log-likelihood function and convergence status when the distribution is fitted to some data}
\usage{
deepg(x, spec, lambda = 1, a = 1, log = FALSE, ...)
peepg(x, spec, lambda = 1, a = 1, log.p = FALSE, lower.tail = TRUE, ...)
qeepg(p, spec, lambda = 1, a = 1, log.p = FALSE, lower.tail = TRUE, ...)
reepg(n, spec, lambda = 1, a = 1, ...)
meepg(g, data, starts, method = "BFGS")
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of probabilities at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{lambda}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the shape parameter, must be positive, the default is 1}
  \item{spec}{a character string specifying the distribution of G and g (for example, "norm" if G and g correspond to the standard normal).}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
  \item{...}{other parameters}
  \item{g}{same as spec but must be one of chisquare ("chisq"), exponential ("exp"), F ("f"), gamma ("gamma"), lognormal ("lognormal"), Weibull ("weibull"), Burr XII ("burrxii"), Chen ("chen"), Frechet ("frechet"), Gompertz ("gompertz"), linear failure rate ("lfr"), log-logistic ("log-logistic"), Lomax ("lomax") and Rayleigh ("rayleigh").  Each of these distributions has one parameter (\code{r}) or two parameters (\code{r}, \code{s}),  for details  including the density function and parameter specifications see Nadarajah and Rocha (2014)}
  \item{data}{a vector of data values for which the distribution is to be fitted}
  \item{starts}{initial values of \code{(lambda, a, r)} if \code{g} has one parameter or initial values of \code{(lambda, a, r, s)} if \code{g} has two parameters}
  \item{method}{the method for optimizing the log likelihood function.  It can be one of \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"} or \code{"SANN"}.  The default is \code{"BFGS"}.  The details of these methods can be found in the manual pages for  \code{optim}}}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random numbers generated or an object giving the values of Cramer-von Misses statistic, Anderson Darling statistic, Kolmogorov Smirnov test statistic and p-value, maximum likelihood estimates, Akaike Information Criterion, Consistent Akaikes Information Criterion, Bayesian Information Criterion, Hannan-Quinn information criterion, standard errors of the maximum likelihood estimates, minimum value of the negative log-likelihood function and convergence status.}
\references{S. Nadarajah and R. Rocha, Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software, 69(10), 1-32, doi:10.18637/jss.v069.i10

M. M. Ristic, S. Nadarajah, A new lifetime distribution, Journal of Statistical Computation and Simulation, doi: 10.1080/00949655.2012.697163}
\author{Saralees Nadarajah, Ricardo Rocha}
\examples{x=runif(10,min=0,max=1)
deepg(x,"exp",lambda=1,a=1)
peepg(x,"exp",lambda=1,a=1)
qeepg(x,"exp",lambda=1,a=1)
reepg(10,"exp",lambda=1,a=1)
meepg("exp",rexp(100),starts=c(1,1,1),method="BFGS")}
\keyword{Exponentiated exponential Poisson G distribution}
