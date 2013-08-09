

#THE GAMMA-UNIFORM DISTRIBUTION AND ITS APPLICATIONS
#H. Torabi and N. H. Montazeri
#Kybernetika, vol. 48 (2012), 16-30


dgammag<-function(x, spec, a=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-f(x, ...)*(1-F(x, ...))**(-2)*dgamma(F(x, ...)/(1-F(x, ...)),shape=a)
        pdf[log==TRUE]<-fL(x, ...)-2*log(1-F(x, ...))+dgamma(F(x, ...)/(1-F(x, ...)),shape=a,log=TRUE)
        return(pdf)

}



pgammag<-function(x, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pgamma(F(x, ...)/(1-F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pgamma(F(x, ...)/(1-F(x, ...)),shape=a,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-1-pgamma(F(x, ...)/(1-F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(1-pgamma(F(x, ...)/(1-F(x, ...)),shape=a))
        return(cdf)

}



qgammag<-function(p, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(qgamma(p,shape=a)/(1+qgamma(p,shape=a)), ...)
	return(qf)

}


rgammag<-function(n, spec, a=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qgammag(u, spec, a=a, ...)
	return(sf)
}



#A new method for generating families of continuous distributions
#A. Alzaatreh, C. Lee, F. Famoye
#METRON (2013) 71 63-79



dbetaexpg<-function(x, spec, lambda=1, a=1, b=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-lambda*f(x, ...)*(1-F(x, ...))**(b-1)*dbeta(1-(1-F(x, ...))**lambda,shape1=a,shape2=b)
        pdf[log==TRUE]<-log(lambda)+fL(x, ...)+(b-1)*log(1-F(x, ...))+dbeta(1-(1-F(x, ...))**lambda,shape1=a,shape2=b,log=TRUE)
        return(pdf)

}


pbetaexpg<-function(x, spec, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-1-pbeta((1-F(x, ...))**lambda,shape1=lambda*(b-1)+1,shape2=a)
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-pbeta((1-F(x, ...))**lambda,shape1=lambda*(b-1)+1,shape2=a))
        cdf[log.p==FALSE&lower.tail==FALSE]<-pbeta((1-F(x, ...))**lambda,shape1=lambda*(b-1)+1,shape2=a)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pbeta((1-F(x, ...))**lambda,shape1=lambda*(b-1)+1,shape2=a,log.p=TRUE)
	return(cdf)

}


qbetaexpg<-function(p, spec, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(1-(qbeta(1-p,shape1=lambda*(b-1)+1,shape2=a))**(1/lambda), ...)
	return(qf)

}


rbetaexpg<-function(n, spec, lambda=1, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qbetaexpg(u, spec, lambda=lambda, a=a, b=b, ...)
	return(sf)
}

#A new method for generating families of continuous distributions
#A. Alzaatreh, C. Lee, F. Famoye
#METRON (2013) 71 63-79


dweibullg<-function(x, spec, beta=1, c=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-c*beta**(-c)*(f(x, ...)/(1-F(x, ...)))*(-log(1-F(x, ...)))**(c-1)*exp(-beta**(-c)*(-log(1-F(x, ...)))**c)
        pdf[log==TRUE]<-log(c)-c*log(beta)+fL(x, ...)-log(1-F(x, ...))+(c-1)*log(-log(1-F(x, ...)))-beta**(-c)*(-log(1-F(x, ...)))**c
        return(pdf)

}


pweibullg<-function(x, spec, beta=1, c=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-1-exp(-beta**(-c)*(-log(1-F(x, ...)))**c)
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-exp(-beta**(-c)*(-log(1-F(x, ...)))**c))
        cdf[log.p==FALSE&lower.tail==FALSE]<-exp(-beta**(-c)*(-log(1-F(x, ...)))**c)
        cdf[log.p==TRUE&lower.tail==FALSE]<--beta**(-c)*(-log(1-F(x, ...)))**c
        return(cdf)

}


qweibullg<-function(p, spec, beta=1, c=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(1-exp(-beta*(-log(1-p))**(1/c)), ...)
	return(qf)

}


rweibullg<-function(n, spec, beta=1, c=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qweibullg(u, spec, beta=beta, c=c, ...)
	return(sf)
}




#The geometric exponential Poisson distribution
#S. Nadarajah, V. G. Cancho, E. M. M. Ortega
#Stat Methods Appl (2013) 22 355-380



dgepg<-function(x, spec, theta=1, eta=0.5, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-theta*(1-eta)*(1-exp(-theta))*f(x, ...)*exp(-theta+theta*F(x, ...))/(1-exp(-theta)-eta+eta*exp(-theta+theta*F(x, ...)))**2
        pdf[log==TRUE]<-log(theta)+log(1-eta)+log(1-exp(-theta))+fL(x, ...)-theta+theta*F(x, ...)-2*log(1-exp(-theta)-eta+eta*exp(-theta+theta*F(x, ...)))
	return(pdf)

}


pgepg<-function(x, spec, theta=1, eta=0.5, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-(exp(-theta+theta*F(x, ...))-exp(-theta))/(1-exp(-theta)-eta+eta*exp(-theta+theta*F(x, ...)))
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(exp(-theta+theta*F(x, ...))-exp(-theta))-log(1-exp(-theta)-eta+eta*exp(-theta+theta*F(x, ...)))
        cdf[log.p==FALSE&lower.tail==FALSE]<-(1-eta)*(1-exp(-theta+theta*F(x, ...)))/(1-exp(-theta)-eta+eta*exp(-theta+theta*F(x, ...)))
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(1-eta)+log(1-exp(-theta+theta*F(x, ...)))-log(1-exp(-theta)-eta+eta*exp(-theta+theta*F(x, ...)))
	return(cdf)

}


qgepg<-function(p, spec, theta=1, eta=0.5, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq((1/theta)*log(((1-exp(-theta)-eta)*p+exp(-theta))/((1-eta*p)*exp(-theta))), ...)
	return(qf)

}


rgepg<-function(n, spec, theta=1, eta=0.5, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qgepg(u, spec, theta=theta, eta=eta, ...)
	return(sf)
}



#A new lifetime distribution
#M. M. Ristic and S. Nadarajah
#Journal of Statistical Computation and Simulation, doi: 10.1080/00949655.2012.697163



deepg<-function(x, spec, lambda=1, a=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	pdf<-x
        pdf[log==FALSE]<-a*lambda*f(x, ...)*(F(x, ...))**(a-1)*exp(-lambda*(F(x, ...))**a)/(1-exp(-lambda))
        pdf[log==TRUE]<-log(lambda)+log(a)+fL(x, ...)+(a-1)*FL(x, ...)-lambda*(F(x, ...))**a-log(1-exp(-lambda))
	return(pdf)

}


peepg<-function(x, spec, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-(1-exp(-lambda*(F(x, ...))**a))/(1-exp(-lambda))
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-exp(-lambda*(F(x, ...))**a))-log(1-exp(-lambda))
        cdf[log.p==FALSE&lower.tail==FALSE]<-(exp(-lambda*(F(x, ...))**a)-exp(-lambda))/(1-exp(-lambda))
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(exp(-lambda*(F(x, ...))**a)-exp(-lambda))-log(1-exp(-lambda))
	return(cdf)

}


qeepg<-function(p, spec, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq((-(1/lambda)*log(1-p*(1-exp(-lambda))))**(1/a), ...)
	return(qf)

}


reepg<-function(n, spec, lambda=1, a=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qeepg(u, spec, lambda=lambda, a=a, ...)
	return(sf)
}


#Truncated-exponential skew-symmetric distributions
#S. Nadarajah, V. Nassiri and A. Mohammadpour
#Statistics, to appear


dtessg<-function(x, spec, lambda=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-lambda*f(x, ...)*exp(-lambda*F(x, ...))/(1-exp(-lambda))
        pdf[log==TRUE]<-log(lambda)+fL(x, ...)-lambda*F(x, ...)-log(1-exp(-lambda))
	return(pdf)

}


ptessg<-function(x, spec, lambda=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-(1-exp(-lambda*F(x, ...)))/(1-exp(-lambda))
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-exp(-lambda*F(x, ...)))-log(1-exp(-lambda))
        cdf[log.p==FALSE&lower.tail==FALSE]<-(exp(-lambda*F(x, ...))-exp(-lambda))/(1-exp(-lambda))
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(exp(-lambda*F(x, ...))-exp(-lambda))-log(1-exp(-lambda))
	return(cdf)

}


qtessg<-function(p, spec, lambda=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(-(1/lambda)*log(1-p*(1-exp(-lambda))), ...)
	return(qf)

}


rtessg<-function(n, spec, lambda=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qtessg(u, spec, lambda=lambda, ...)
	return(sf)
}


#Marshall, A. W. and Olkin, I. (1997).
#A new method for adding a parameter to a family of distributions with application to the exponential and Weibull families.
#Biometrika, 84, 641-652.




dmog<-function(x, spec, beta=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-beta*f(x, ...)/(1-(1-beta)*F(x, ...))**2
        pdf[log==TRUE]<-log(beta)+fL(x, ...)-2*log(1-(1-beta)*F(x, ...))
	return(pdf)

}

pmog<-function(x, spec, beta=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-beta*F(x, ...)/(1-(1-beta)*F(x, ...))
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(beta)+FL(x, ...)-log(1-(1-beta)*F(x, ...))
        cdf[log.p==FALSE&lower.tail==FALSE]<-(1-F(x, ...))/(1-(1-beta)*F(x, ...))
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(1-F(x, ...))-log(1-(1-beta)*F(x, ...))
	return(cdf)

}


qmog<-function(p, spec, beta=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(p/(beta+(1-beta)*p))
	return(qf)

}

rmog<-function(n, spec, beta=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qmog(u, spec, beta=beta, ...)
	return(sf)
}


#Gupta, R. C., Gupta, P. L. and Gupta, R. D.  (1998).
#Modeling failure time data by Lehman alternatives.
#Communications in Statistics---Theory and Methods, 27, 887-904.



dexpg<-function(x, spec, a=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	pdf<-x
        pdf[log==FALSE]<-a*f(x, ...)*(F(x, ...))**(a-1)
        pdf[log==TRUE]<-log(a)+fL(x, ...)+(a-1)*FL(x, ...)
	return(pdf)

}


pexpg<-function(x, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-(F(x, ...))**a
        cdf[log.p==TRUE&lower.tail==TRUE]<-a*FL(x, ...)
        cdf[log.p==FALSE&lower.tail==FALSE]<-1-(F(x, ...))**a
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(1-(F(x, ...))**a)
	return(cdf)

}


qexpg<-function(p, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(p**(1/a), ...)
	return(qf)

}


rexpg<-function(n, spec, a=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qexpg(u, spec, a=a, ...)
	return(sf)
}


#Lemonte, A. J., Barreto-Souza, W. and Cordeiro, G. M. (2013).
#The exponentiated Kumaraswamy distribution and its log-transform.
#Brazilian Journal of Probability and Statistics, 27, 31-53. 


dexpkumg<-function(x, spec, a=1, b=1, c=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	pdf<-x
        pdf[log==FALSE]<-a*b*c*f(x, ...)*(F(x, ...))**(a-1)*(1-(F(x, ...))**a)**(b-1)*(1-(1-(F(x, ...))**a)**b)**(c-1)
        pdf[log==TRUE]<-log(a)+log(b)+log(c)+fL(x, ...)+(a-1)*FL(x, ...)+(b-1)*log(1-(F(x, ...))**a)+(c-1)*log(1-(1-(F(x, ...))**a)**b)
	return(pdf)

}


pexpkumg<-function(x, spec, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-(1-(1-(F(x, ...))**a)**b)**c
        cdf[log.p==TRUE&lower.tail==TRUE]<-c*log(1-(1-(F(x, ...))**a)**b)
        cdf[log.p==FALSE&lower.tail==FALSE]<-1-(1-(1-(F(x, ...))**a)**b)**c
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(1-(1-(1-(F(x, ...))**a)**b)**c)
	return(cdf)

}


qexpkumg<-function(p, spec, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq((1-(1-p**(1/c))**(1/b))**(1/a), ...)
	return(qf)

}


rexpkumg<-function(n, spec, a=1, b=1, c=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qexpkumg(u, spec, a=a, b=b, c=c, ...)
	return(sf)
}


#Cordeiro, Gauss M.; Ortega, Edwin M. M.; da Cunha, Daniel C. C.
#The exponentiated generalized class of distributions.
#Journal of Data Science 11 (2013), 1–27.


deg<-function(x, spec, a=1, b=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-a*b*f(x, ...)*(1-F(x, ...))**(a-1)*(1-(1-F(x, ...))**a)**(b-1)
        pdf[log==TRUE]<-log(a)+log(b)+fL(x, ...)+(a-1)*log(1-F(x, ...))+(b-1)*log(1-(1-F(x, ...))**a)
	return(pdf)

}


peg<-function(x, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-(1-(1-F(x, ...))**a)**b
        cdf[log.p==TRUE&lower.tail==TRUE]<-b*log(1-(1-F(x, ...))**a)
        cdf[log.p==FALSE&lower.tail==FALSE]<-1-(1-(1-F(x, ...))**a)**b
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(1-(1-(1-F(x, ...))**a)**b)
	return(cdf)

}


qeg<-function(p, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(1-(1-p**(1/b))**(1/a), ...)
	return(qf)

}

reg<-function(n, spec, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qeg(u, spec, a=a, b=b, ...)
	return(sf)
}


#Cordeiro, G. M. and Castro, M. (2011).
#A new family of generalized distributions.
#Journal of Statistical Computation and Simulation, 81, 883-898.



dkumg<-function(x, spec, a=1, b=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	pdf<-x
        pdf[log==FALSE]<-a*b*f(x, ...)*(F(x, ...))**(a-1)*(1-F(x, ...))**(b-1)
        pdf[log==TRUE]<-log(a)+log(b)+fL(x, ...)+(a-1)*FL(x, ...)+(b-1)*log(1-F(x, ...))
	return(pdf)

}



pkumg<-function(x, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-1-(1-(F(x, ...))**a)**b
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-(1-(F(x, ...))**a)**b)
        cdf[log.p==FALSE&lower.tail==FALSE]<-(1-(F(x, ...))**a)**b
        cdf[log.p==TRUE&lower.tail==FALSE]<-b*log(1-(F(x, ...))**a)
	return(cdf)

}


qkumg<-function(p, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq((1-(1-p)**(1/b))**(1/a), ...)
	return(qf)

}


rkumg<-function(n, spec, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qkumg(u, spec, a=a, b=b, ...)
	return(sf)
}


#Cordeiro, G. M., Ortega, E. M. M. and Silva, G. (2012a).
#The beta extended Weibull family.
#Journal of Probability and Statistical Science, 10, 15-40.



dbeg<-function(x, spec, alpha=1, a=1, b=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-f(x, ...)*exp(-alpha*F(x, ...))*dbeta(exp(-alpha*F(x, ...)),shape1=a,shape2=b)
        pdf[log==TRUE]<-fL(x, ...)-alpha*F(x, ...)+dbeta(exp(-alpha*F(x, ...)),shape1=a,shape2=b,log=TRUE)
	return(pdf)

}


pbeg<-function(x, spec, alpha=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pbeta(1-exp(-alpha*F(x, ...)),shape1=a,shape2=b)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pbeta(1-exp(-alpha*F(x, ...)),shape1=a,shape2=b,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-pbeta(1-exp(-alpha*F(x, ...)),shape1=a,shape2=b,lower.tail==FALSE)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pbeta(1-exp(-alpha*F(x, ...)),shape1=a,shape2=b,log.p=TRUE,lower.tail==FALSE)
	return(cdf)

}


qbeg<-function(p, spec, alpha=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(-(1/alpha)*log(1-qbeta(p,shape1=a,shape2=b)), ...)
	return(qf)

}


rbeg<-function(n, spec, alpha=1, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qbeg(u, spec, alpha=alpha, a=a, b=b, ...)
	return(sf)
}


#Alexander, C., Cordeiro, G. M., Ortega, E. M. M. and Sarabia, J. M. (2012).
#Generalized beta-generated distributions.
#Computational Statistics and Data Analysis, 56, 1880-1897.



dgbg<-function(x, spec, a=1, b=1, c=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	pdf<-x
        pdf[log==FALSE]<-c*f(x, ...)*(F(x, ...))**(c-1)*dbeta((F(x, ...))**c,shape1=a,shape2=b)
        pdf[log==TRUE]<-log(c)+fL(x, ...)+(c-1)*FL(x, ...)+dbeta((F(x, ...))**c,shape1=a,shape2=b,log=TRUE)
	return(pdf)

}


pgbg<-function(x, spec, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pbeta((F(x, ...))**c,shape1=a,shape2=b)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pbeta((F(x, ...))**c,shape1=a,shape2=b,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-pbeta((F(x, ...))**c,shape1=a,shape2=b,lower.tail==FALSE)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pbeta((F(x, ...))**c,shape1=a,shape2=b,log.p=TRUE,lower.tail==FALSE)
	return(cdf)

}


qgbg<-function(p, spec, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq((qbeta(p,shape1=a,shape2=b))**(1/c), ...)
	return(qf)

}


rgbg<-function(n, spec, a=1, b=1, c=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qgbg(u, spec, a=a, b=b, c=c, ...)
	return(sf)
}


#Modified beta distributions
#S. Nadarajah, M. Teimouri, S. H. Shih
#Sankhya, to appear


dmbetag<-function(x, spec, beta=1, a=1, b=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-beta**a*f(x, ...)*dbeta(F(x, ...),shape1=a,shape2=b)/(1-(1-beta)*F(x, ...))**(a+b)
        pdf[log==TRUE]<-a*log(beta)+fL(x, ...)+dbeta(F(x, ...),shape1=a,shape2=b,log=TRUE)-(a+b)*log(1-(1-beta)*F(x, ...))
	return(pdf)

}




pmbetag<-function(x, spec, beta=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pbeta(beta*F(x, ...)/(1-(1-beta)*F(x, ...)),shape1=a,shape2=b)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pbeta(beta*F(x, ...)/(1-(1-beta)*F(x, ...)),shape1=a,shape2=b,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-pbeta(beta*F(x, ...)/(1-(1-beta)*F(x, ...)),shape1=a,shape2=b,lower.tail==FALSE)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pbeta(beta*F(x, ...)/(1-(1-beta)*F(x, ...)),shape1=a,shape2=b,log.p=TRUE,lower.tail==FALSE)
	return(cdf)

}



qmbetag<-function(p, spec, beta=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(qbeta(p,shape1=a,shape2=b)/(beta-(beta-1)*qbeta(p,shape1=a,shape2=b)), ...)
	return(qf)

}


rmbetag<-function(n, spec, beta=1, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qmbetag(u, spec, beta=beta, a=a, b=b, ...)
	return(sf)
}


#Eugene, N., Lee, C. and Famoye, F. (2002).
#Beta-normal distribution and its applications.
#Communications in Statistics---Theory and Methods, 31, 497-512.



dbetag<-function(x, spec, a=1, b=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-f(x, ...)*dbeta(F(x, ...),shape1=a,shape2=b)
        pdf[log==TRUE]<-fL(x, ...)+dbeta(F(x, ...),shape1=a,shape2=b,log=TRUE)
	return(pdf)

}


pbetag<-function(x, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pbeta(F(x, ...),shape1=a,shape2=b)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pbeta(F(x, ...),shape1=a,shape2=b,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-pbeta(F(x, ...),shape1=a,shape2=b,lower.tail==FALSE)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pbeta(F(x, ...),shape1=a,shape2=b,log.p=TRUE,lower.tail==FALSE)
	return(cdf)

}

qbetag<-function(p, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(qbeta(p,shape1=a,shape2=b), ...)
	return(qf)

}


rbetag<-function(n, spec, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qbetag(u, spec, a=a, b=b, ...)
	return(sf)
}


#M. Amini, S. M. T. K. MirMostafaee, J. Ahmadi
#Log-gamma-generated families of distributions
#Statistics, doi: 10.1080/02331888.2012.748775



dloggammag1<-function(x, spec, a=1, b=1, log=FALSE, ...)
{

        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-b**a*f(x, ...)*(-log(1-F(x, ...)))**(a-1)*(1-F(x, ...))**(b-1)/gamma(a)
        pdf[log==TRUE]<-a*log(b)+fL(x, ...)+(a-1)*log(-log(1-F(x, ...)))+(b-1)*log(1-F(x, ...))-lgamma(a)
	return(pdf)

}


ploggammag1<-function(x, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{

        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-1-pgamma(-b*log(1-F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-pgamma(-b*log(1-F(x, ...)),shape=a))
        cdf[log.p==FALSE&lower.tail==FALSE]<-pgamma(-b*log(1-F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pgamma(-b*log(1-F(x, ...)),shape=a,log.p=TRUE)
	return(cdf)

}


qloggammag1<-function(p, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(1-exp(-(1/b)*qgamma(1-p,shape=a)), ...)
	return(qf)

}


rloggammag1<-function(n, spec, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qloggammag1(u, spec, a=a, b=b, ...)
	return(sf)
}

#M. Amini, S. M. T. K. MirMostafaee, J. Ahmadi
#Log-gamma-generated families of distributions
#Statistics, doi: 10.1080/02331888.2012.748775




dloggammag2<-function(x, spec, a=1, b=1, log=FALSE, ...)
{

        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	pdf<-x
        pdf[log==FALSE]<-b**a*f(x, ...)*(-FL(x, ...))**(a-1)*(F(x, ...))**(b-1)/gamma(a)
        pdf[log==TRUE]<-a*log(b)+fL(x, ...)+(a-1)*log(-FL(x, ...))+(b-1)*FL(x, ...)-lgamma(a)
	return(pdf)

}


ploggammag2<-function(x, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{

        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pgamma(-b*FL(x, ...),shape=a)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pgamma(-b*FL(x, ...),shape=a,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-1-pgamma(-b*FL(x, ...),shape=a)
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(1-pgamma(-b*FL(x, ...),shape=a))
	return(cdf)

}


qloggammag2<-function(p, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(exp(-(1/b)*qgamma(p,shape=a)), ...)
	return(qf)

}


rloggammag2<-function(n, spec, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qloggammag2(u, spec, a=a, b=b, ...)
	return(sf)
}


#Zografos, K. and Balakrishnan, N. (2009).
#On families of beta- and generalized gamma-generated distributions and associated inference.
#Statistical Methodology, 6, 344-362.




dgammag1<-function(x, spec, a=1, log=FALSE, ...)
{

        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-f(x, ...)*(-log(1-F(x, ...)))**(a-1)/gamma(a)
        pdf[log==TRUE]<-fL(x, ...)+(a-1)*log(-log(1-F(x, ...)))-lgamma(a)
	return(pdf)

}


pgammag1<-function(x, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{

        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pgamma(-log(1-F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pgamma(-log(1-F(x, ...)),shape=a,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-pgamma(-log(1-F(x, ...)),shape=a,lower.tail==FALSE)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pgamma(-log(1-F(x, ...)),shape=a,log.p=TRUE,lower.tail==FALSE)
	return(cdf)

}

qgammag1<-function(p, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(1-exp(-qgamma(p,shape=a)), ...)
	return(qf)

}


rgammag1<-function(n, spec, a=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qgammag1(u, spec, a=a, ...)
	return(sf)
}


#Ristic, M. M. and Balakrishnan, N. (2012).
#The gamma exponentiated exponential distribution.
#Journal of Statistical Computation and Simulation, 82, 1191-1206.


dgammag2<-function(x, spec, a=1, log=FALSE, ...)
{

        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-f(x, ...)*(-log(F(x, ...)))**(a-1)/gamma(a)
        pdf[log==TRUE]<-fL(x, ...)+(a-1)*log(-log(F(x, ...)))-lgamma(a)
	return(pdf)

}

pgammag2<-function(x, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{

        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-1-pgamma(-log(F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-pgamma(-log(F(x, ...)),shape=a))
        cdf[log.p==FALSE&lower.tail==FALSE]<-pgamma(-log(F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pgamma(-log(F(x, ...)),shape=a,log.p=TRUE)
	return(cdf)

}

qgammag2<-function(p, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-p
        qf<-Fq(exp(-qgamma(1-p,shape=a)), ...)
	return(qf)

}

rgammag2<-function(n, spec, a=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qgammag2(u, spec, a=a, ...)
	return(sf)
}
