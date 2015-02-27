\name{tessg}
\alias{dtessg}
\alias{ptessg}
\alias{qtessg}
\alias{rtessg}
\alias{mtessg}
\title{Truncated-Exponential Skew-Symmetric G Distribution}
\description{Computes the pdf, cdf, quantile and random numbers of the truncated-exponential skew-symmetric G distribution due to Nadarajah et al. (2013) specified by the pdf
\deqn{f (x) = \frac {\lambda}{1 - \exp (-\lambda)} g(x) \exp \left\{ -\lambda G(x) \right\}}
for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf, and \eqn{-\infty < \lambda < \infty}, the skewness parameter.  Also computes the Cramer-von Misses statistic, Anderson Darling statistic, Kolmogorov Smirnov test statistic and p-value, maximum likelihood estimates, Akaike Information Criterion, Consistent Akaikes Information Criterion, Bayesian Information Criterion, Hannan-Quinn information criterion, standard errors of the maximum likelihood estimates, minimum value of the negative log-likelihood function and convergence status when the distribution is fitted to some data}
\usage{
dtessg(x, spec, lambda = 1, log = FALSE, ...)
ptessg(x, spec, lambda = 1, log.p = FALSE, lower.tail = TRUE, ...)
qtessg(p, spec, lambda = 1, log.p = FALSE, lower.tail = TRUE, ...)
rtessg(n, spec, lambda = 1, ...)
mtessg(g, data, starts, method = "BFGS")
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of probabilities at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{lambda}{the value of skewness parameter, can be any real value, the default is 1}
  \item{spec}{a character string specifying the distribution of G and g (for example, "norm" if G and g correspond to the standard normal).}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
  \item{...}{other parameters}
  \item{g}{same as spec but must be one of chisquare ("chisq"), exponential ("exp"), F ("f"), gamma ("gamma"), lognormal ("lognormal"), Weibull ("weibull"), Burr XII ("burrxii"), Chen ("chen"), Frechet ("frechet"), Gompertz ("gompertz"), linear failure rate ("lfr"), log-logistic ("log-logistic"), Lomax ("lomax") and Rayleigh ("rayleigh").  Each of these distributions has one parameter (\code{r}) or two parameters (\code{r}, \code{s}),  for details  including the density function and parameter specifications see Nadarajah and Rocha (2014)}
  \item{data}{a vector of data values for which the distribution is to be fitted}
  \item{starts}{initial values of \code{(lambda, r)} if \code{g} has one parameter or initial values of \code{(lambda, r, s)} if \code{g} has two parameters}
  \item{method}{the method for optimizing the log likelihood function.  It can be one of \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"} or \code{"SANN"}.  The default is \code{"BFGS"}.  The details of these methods can be found in the manual pages for  \code{optim}}}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random numbers generated or an object giving the values of Cramer-von Misses statistic, Anderson Darling statistic, Kolmogorov Smirnov test statistic and p-value, maximum likelihood estimates, Akaike Information Criterion, Consistent Akaikes Information Criterion, Bayesian Information Criterion, Hannan-Quinn information criterion, standard errors of the maximum likelihood estimates, minimum value of the negative log-likelihood function and convergence status.}
\references{S. Nadarajah and R. Rocha, Newdistns: An R Package for new families of distributions, submitted 

S. Nadarajah, V. Nassiri, A. Mohammadpour, Truncated-exponential skew-symmetric distributions, Statistics, to appear}
\author{Saralees Nadarajah, Ricardo Rocha}
\examples{x=runif(10,min=0,max=1)
dtessg(x,"exp",lambda=1)
ptessg(x,"exp",lambda=1)
qtessg(x,"exp",lambda=1)
rtessg(10,"exp",lambda=1)
mtessg("exp",rexp(100),starts=c(1,1),method="BFGS")}
\keyword{Truncated-exponential skew-symmetric G distribution}
