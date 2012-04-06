\name{condhparetomixt}
\alias{condhparetomixt}
\alias{pcondhparetomixt}
\alias{pcondhparetomixt.dirac}
\alias{dcondhparetomixt}
\title{
The conditional hybrid Pareto mixture distribution
}
\description{
Distribution function for the conditional hybrid Pareto mixture with and
without discrete component at zero and density function for the
conditional hybrid Pareto mixture without discrete component.
}
\usage{
pcondhparetomixt(params, m, y, trunc = TRUE)
pcondhparetomixt.dirac(params,m,y)
dcondhparetomixt(params,m,y,log=FALSE,trunc=TRUE)
}
\arguments{
  \item{params}{\code{m} x 4 x n matrix of mixture parameters where n
    is the number of examples for \code{condhparetomixt} and
    (4\code{m}+1) x n matrix for \code{condhparetomixt.dirac}
}
  \item{m}{
    Number of mixture components.
}
  \item{y}{
    Vector of n dependent variables.
  }
  \item{log}{logical, if TRUE, probabilities p are
    given as log(p).}
  \item{trunc}{logical, if TRUE (default), the hybrid Pareto density is
    truncated below zero. A mixture with a Dirac component at zero is
    always truncated below zero.
    }
}
\details{
\code{params} can be computed from the forward functions on the
explanatory variables x of dimension d x n associated with \code{y} :
\code{condhparetomixt.fwd} and \code{condhparetomixt.dirac.fwd}
}
\value{
  Distribution function evaluated at n points for
  \code{pcondhparetomixt} and \code{pcondhparetomixt.dirac} and density
  function for \code{dcondhparetomixt}.
}
\references{
  Carreau, J. and Bengio, Y. (2009), A Hybrid {P}areto Mixture for
  Conditional Asymmetric Fat-Tailed Distributions, 20, IEEE Transactions
  on Neural Networks
}
\author{
Julie Carreau
}
\seealso{
\code{\link{condmixt.nll}}, \code{\link{condmixt.fwd}}
}
\examples{
# generate train data with a mass at zero
ntrain <- 200
xtrain <- runif(ntrain,0,2*pi)
alpha.train <- sin(2*xtrain)/2+1/2
data.train <- rep(0,ntrain)
for (i in 1:ntrain){
  if (sample(c(0,1),1,prob=c(1-alpha.train[i],alpha.train[i]))){
   # rainy day, sample from a Frechet
    data.train[i] <-rfrechet(1,loc=3*sin(2*xtrain[i])+4,scale=1/(1+exp(-(xtrain[i]-1))),
                    shape=(1+exp(-(xtrain[i]/5-2))))
  }
}

plot(xtrain,data.train,pch=20)


h <- 4 # number of hidden units
m <- 2 # number of components

# initialize a conditional mixture with hybrid Pareto components and a dirac at zero
thetainit <- condhparetomixt.dirac.init(1,h,m,data.train)


# compute mixture parameters on test data
params.mixt <- condhparetomixt.dirac.fwd(thetainit,h,m,t(xtrain))



cdf <- pcondhparetomixt.dirac(params.mixt,m,data.train) # compute CDF on test data

}
