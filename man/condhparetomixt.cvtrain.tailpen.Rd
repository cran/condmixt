\name{condhparetomixt.cvtrain.tailpen}
\alias{condhparetomixt.cvtrain.tailpen}
\title{
Cross-validation of the conditinal mixture with hybrid Pareto components
with a tail penalty added to the negative log-likelihood for training.
}
\description{
K-fold cross-validation is performed in order to select hyper-parameters
of the conditional mixture with hybrid Pareto components. A tail penalty
is added to the negative log-likelihood in order to guide the tail index
parameters estimation. Performance is evaluated by computing the
negative log-likelihood, without penalty.
}
\usage{
condhparetomixt.cvtrain.tailpen(x, y, hp, nfold = 5, nstart = 1, ...)
}
\arguments{
  \item{x}{
    Matrix of explanatory (independent) variables of dimension d x n, d is
    the number of variables and n is the number of examples (patterns)
}
\item{y}{
  Vector of n dependent variables
}
  \item{hp}{Matrix nhp x 7 whose rows represent a set of
    values for the following hyper-parameters : number of hidden unit,
    number of component, lambda, w, beta, mu and sigma. The last five
    hyper-parameters control the tail penalty, see in
  the Details section below.
}
  \item{nfold}{
Number of folds for the cross-validation, default is 5.
}
  \item{nstart}{
Number of optimization re-starts for each training, default is
one. Optimization is re-started from different initial values in order
to avoir local minima.
}
  \item{\dots}{
Extra arguments passed to \code{nlm}.
}
}
\details{
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
\value{
Returns a vector of negative log-likelihood values evaluated
out-of-sample by cross-validation. Elements in the vector correspond to
each set of hyper-parameters evaluated.
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
\code{\link{condmixt.foldtrain}},\code{\link{condmixt.train}},\code{\link{condmixt.nll}}, \code{\link{condmixt.init}}
}
\examples{
n <- 200
x <- runif(n) # x is a random uniform variate
 # y depends on x through the parameters of the Frechet distribution
y <- rfrechet(n,loc = 3*x+1,scale = 0.5*x+0.001,shape=x+1)

plot(x,y,pch=22)
# 0.99 quantile of the generative distribution
qgen <- qfrechet(0.99,loc = 3*x+1,scale = 0.5*x+0.001,shape=x+1)
points(x,qgen,pch="*",col="orange")

# create a matrix with sets of values for the number of hidden units and
# the number of components
hp <- matrix(c(2,4,2,2),nrow=2,ncol=2)

# keep tail penalty parameters constant
hp <- cbind(hp, rep(10,2),rep(0.5,2),rep(20,2),rep(0.2,2),rep(0.5,2))

condhparetomixt.cvtrain.tailpen(t(x), y, hp, nfold = 2, nstart = 2, iterlim=100)
}
