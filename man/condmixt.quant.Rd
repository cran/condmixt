\name{condmixt.quant}
\alias{condmixt.quant}
\alias{condhparetomixt.quant}
\alias{condhparetomixt.dirac.quant}
\alias{condhparetomixt.dirac.condquant}
\alias{condgaussmixt.quant}
\alias{condgaussmixt.dirac.quant}
\alias{condgaussmixt.dirac.condquant}
\alias{condlognormixt.quant}
\alias{condlognormixt.dirac.quant}
\alias{condlognormixt.dirac.condquant}
\alias{condbergamixt.quant}
\title{
Quantile computation for conditional mixtures.
}
\description{
Quantile computation for conditional mixtures requires to solve
numerically F(y)=p where F is the distribution function of the
conditional mixture and p is a probability level.
}
\usage{
condhparetomixt.quant(theta, h, m, x, p, a, b, trunc = TRUE)
condhparetomixt.dirac.quant(theta,h,m,x,p,a,b)
condhparetomixt.dirac.condquant(theta,h,m,x,p,a,b)
condgaussmixt.quant(theta,h,m,x,p,a,b,trunc=TRUE)
condgaussmixt.dirac.quant(theta,h,m,x,p,a,b)
condgaussmixt.dirac.condquant(theta,h,m,x,p,a,b)
condlognormixt.quant(theta,h,m,x,p,a,b)
condlognormixt.dirac.quant(theta,h,m,x,p,a,b)
condlognormixt.dirac.condquant(theta,h,m,x,p,a,b)
condbergamixt.quant(theta,h,x,p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
    Vector of neural network parameters
}
  \item{h}{
    Number of hidden units
}
  \item{m}{
    Number of components
}
  \item{x}{
    Matrix of explanatory (independent) variables of dimension d x n, d is
    the number of variables and n is the number of examples (patterns)
}
  \item{p}{
    Probability level in [0,1]
}
  \item{a}{
    Approximate lower bound on quantile value.
}
  \item{b}{
    Approximate upper bound on quantile value.
}
  \item{trunc}{
    Logical variable, if true, density is truncated below zero and
    re-weighted to make sure it integrates to one.
}
}
\details{
\code{condhparetomixt} indicates a mixture with hybrid Pareto
components,
\code{condgaussmixt} for Gaussian components,
\code{condlognormixt} for Log-Normal components,
\code{condbergam} for a Bernoulli-Gamma two component mixture,
\code{dirac} indicates that a discrete dirac component is included in
the mixture
\code{condquant} applies for mixtures with a dirac component at zero :
quantiles are computed given that the variable is strictly positive,
that is the quantile is computed for the continuous part of the mixture
only : P(Y <= y | Y >0, X)
}
\value{Computed quantiles are stored in a  matrix whose rows correspond
  to the probability levels and whose columns correspond to the number
  of examples n.
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
\code{\link{condmixt.train}},\code{\link{condmixt.nll}}, \code{\link{condmixt.init}}
}
\examples{
# generate train data
ntrain <- 200
xtrain <- runif(ntrain)
ytrain <- rfrechet(ntrain,loc = 3*xtrain+1,scale =
0.5*xtrain+0.001,shape=xtrain+2)
plot(xtrain,ytrain,pch=22) # plot train data
qgen <- qfrechet(0.99,loc = 3*xtrain+1,scale = 0.5*xtrain+0.001,shape=xtrain+2)
points(xtrain,qgen,pch="*",col="orange")

# generate test data
ntest <- 200
xtest <- runif(ntest)
ytest <- rfrechet(ntest,loc = 3*xtest+1,scale =
0.5*xtest+0.001,shape=xtest+2)

h <- 2 # number of hidden units
m <- 4 # number of components

# train a mixture with hybrid Pareto components
thetaopt <- condhparetomixt.train(h,m,t(xtrain),ytrain, nstart=2,iterlim=100)
qmod <- condhparetomixt.quant(thetaopt,h,m,t(xtest),0.99,0,10,trunc=TRUE)
points(xtest,qmod,pch="o",col="blue")

}
