\name{condmixt.foldtrain}
\alias{condmixt.foldtrain}
\alias{condhparetomixt.foldtrain}
\alias{condhparetomixt.foldtrain.tailpen}
\alias{condhparetomixt.dirac.foldtrain.tailpen}
\alias{condgaussmixt.foldtrain}
\alias{condgaussmixt.dirac.foldtrain}
\alias{condlognormixt.foldtrain}
\alias{condlognormixt.dirac.foldtrain}
\alias{condbergamixt.foldtrain}
\title{
Training of conditional mixtures and evaluation of the negative
log-likelihood on validation data
}
\description{
  This function can be used to parallelize n-fold
  cross-validation. Training is done on a training set and the negative
  log-likelihood is evaluated on validation data. This is repeated for a
  set of hyper-parameter values. Model selection can be performed based
  on the evaluation of each set of hyper-parameters on validation data.
}
\usage{
condhparetomixt.foldtrain(xtrain, ytrain, xtest, ytest, hp, nstart = 1, ...)
condhparetomixt.foldtrain.tailpen(xtrain, ytrain, xtest, ytest, hp, nstart = 1, ...)
condhparetomixt.dirac.foldtrain.tailpen(xtrain, ytrain, xtest, ytest, hp, nstart = 1, ...)
condgaussmixt.foldtrain(xtrain, ytrain, xtest, ytest, hp, nstart = 1, ...)
condgaussmixt.dirac.foldtrain(xtrain, ytrain, xtest, ytest, hp, nstart = 1, ...)
condlognormixt.foldtrain(xtrain, ytrain, xtest, ytest, hp, nstart = 1, ...)
condlognormixt.dirac.foldtrain(xtrain, ytrain, xtest, ytest, hp, nstart = 1, ...)
condbergamixt.foldtrain(xtrain, ytrain, xtest, ytest, hp,nstart = 1, ...)
}
\arguments{
  \item{xtrain}{
Matrix of explanatory (independent) variables of dimension d x ntrain, d is
the number of variables and ntrain is the number of training examples (patterns)
}
  \item{ytrain}{
Vector of ntrain dependent variables
}

  \item{xtest}{
Matrix of explanatory (independent) variables of dimension d x ntest, d is
the number of variables and ntest is the number of test/validation examples (patterns)
}
  \item{ytest}{
Vector of ntest dependent variables
}

\item{hp}{Matrix of hyper-parameters where the columns represent the
  different hyper-parameters and the rows the sets of values for each
  hyper-parameter. For \code{condhparetomixt}, \code{condgaussmixt} and
  \code{lognormixt} together with the version with a dirac component,
  the hyper-parameters are the number of hidden units and the number of
  components. When a tail penalty is included for the hybrid Pareto
  conditional mixture, there are possibly five other hyper-parameters:
  lambda, w, beta, mu and sigma which controls the tail penalty, see in
  the Details section below.

  }
\item{nstart}{Number of minimization re-starts, default is one.}


\item{...}{optional arguments for the optimizer \code{nlm}}
}
\details{
\code{condhparetomixt} indicates a mixture with hybrid Pareto
components,
\code{condgaussmixt} for Gaussian components,
\code{condlognormixt} for Log-Normal components,
\code{condbergam} for a Bernoulli-Gamma two component mixture,
\code{tailpen} indicates that a penalty is added to the log-likelihood
to guide the tail index parameter estimation,
\code{dirac} indicates that a discrete dirac component at zero is included in
the mixture

In order to drive the tail index estimation, a penalty is introduced in
the log-likelihood. The goal of the penalty is to include a priori
information which in our case is that only a few mixture components have
a heavy tail index which should approximate the tail of the underlying
distribution while most other mixture components have a light tail and
aim at modelling the central part of the underlying distribution.

The penalty term is given by the logarithm of the following
two-component mixture, as a function of a tail index parameter xi :

w beta exp(-beta xi) + (1-w) exp(-(xi-mu)^2/(2 sigma^2))/(sqrt(2 pi) sigma)

where the first term is the prior on the light tail component and the
second term is the prior on the heavy tail component.

The extra hyper-parameters for the penalty terms are as follows :
- lambda : penalty parameter which controls the trade-off between
    the penalty and the negative log-likelihood, takes on positive
    values. If zero, no penalty
- w : penalty parameter in [0,1] which is the proportion of
    components with light tails, 1-\code{w} being the proportion of
    components with heavy tails
- beta : positive penalty parameter which indicates the importance
    of the light tail components (it is the parameter of an exponential
    which represents the prior over the light tail components)
- mu : penalty parameter in (0,1) which represents the a priori
    value for the heavy tail index of the underlying distribution
- sigma : positive penalty parameter which controls the spread
  around the a priori value for the heavy tail index of the underlying
  distribution
}
\value{Returns a vector of negative log-likelihood values evaluated on
  the test set corresponding to each set of hyper-parameters.
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
0.5*xtrain+0.001,shape=xtrain+1)
plot(xtrain,ytrain,pch=22) # plot train data
qgen <- qfrechet(0.99,loc = 3*xtrain+1,scale = 0.5*xtrain+0.001,shape=xtrain+1)
points(xtrain,qgen,pch="*",col="orange")

# generate test data
ntest <- 200
xtest <- runif(ntest)
ytest <- rfrechet(ntest,loc = 3*xtest+1,scale =
0.5*xtest+0.001,shape=xtest+1)

# create a matrix with sets of values for the number of hidden units and
# the number of components
hp <- matrix(c(2,4,2,2),nrow=2,ncol=2)

# train and test a mixture with hybrid Pareto components
condhparetomixt.foldtrain(t(xtrain),ytrain,t(xtest),ytest,hp,nstart=2,iterlim=100)
}
