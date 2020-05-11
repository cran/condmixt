\name{condmixt.dirac.negloglike}
\alias{condmixt.dirac.negloglike}
\alias{condhparetomixt.dirac.negloglike}
\alias{condlognormixt.dirac.negloglike}
\alias{condgaussmixt.dirac.negloglike}
\alias{condbergamixt.negloglike}
\title{
Negative log-likelihood for conditional mixture with a discrete
component at zero.
}
\description{
The negative log-likelihood is computed from mixture parameters rather
than from neural network parameters. The mixture parameters are obtained
from the neural network for given explanatory variable values.
}
\usage{
condhparetomixt.dirac.negloglike(params, m, y)
condlognormixt.dirac.negloglike(params,m,y)
condgaussmixt.dirac.negloglike(params,m,y)
condbergamixt.negloglike(params,y)
}
\arguments{
  \item{params}{p x n matrix of mixture parameters where n is the number
    of examples and p = (4\code{m}+1) for \code{condhparetomixt.dirac}, p =
    (3\code{m}+1) for \code{condgaussmixt.dirac} and
    \code{condlognormixt.dirac} and p = 3 for \code{condbergamixt}.
}
  \item{m}{
    Number of components in the mixture.
}
  \item{y}{
Vector of n dependent variables.
}
}
\details{
\code{params} can be computed from the forward functions on the
explanatory variables x of dimension d x n associated with \code{y} :
\code{condhparetomixt.dirac.fwd}, \code{condgaussmixt.dirac.fwd} (which
can be used for conditional mixtures with Log-Normal components) and
\code{condbergamixt.fwd}
}
\value{
Vector of length n corresponding to the negative log-likelihood
evaluated on each example.
}
\references{
Bishop, C. (1995), Neural Networks for Pattern Recognition, Oxford

Carreau, J. and Bengio, Y. (2009), A Hybrid {P}areto Mixture for
Conditional Asymmetric Fat-Tailed Distributions, 20, IEEE Transactions
on Neural Networks
}
\author{
Julie Carreau
}
\seealso{
\code{\link{condmixt.fwd}}
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

# initialize a conditional mixture with hybrid Pareto components and a
# dirac at zero
 thetainit <- condhparetomixt.dirac.init(1,h,m,data.train)

# compute mixture parameters 
params.mixt <- condhparetomixt.dirac.fwd(thetainit,h,m,t(xtrain))

# compute negative log-likelihood 
nll <- condhparetomixt.dirac.negloglike(params.mixt, m, data.train)

}
