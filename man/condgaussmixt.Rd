\name{condgaussmixt}
\alias{condgaussmixt}
\alias{pcondgaussmixt}
\alias{dcondgaussmixt}
\title{
The conditional Gaussian mixture distribution
}
\description{
Distribution function and density function for the conditional Gaussian
mixture without discrete component.
}
\usage{
pcondgaussmixt(params, m, y, trunc = TRUE)
dcondgaussmixt(params,m,y,log=FALSE,trunc=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{\code{m} x 3 x n matrix of mixture parameters where n
    is the number of examples
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
    truncated below zero.
    }
}
\details{
\code{params} can be computed by applying \code{condgaussmixt.fwd} on the
explanatory variables x of dimension d x n associated with \code{y}.
}
\value{
  Distribution function evaluated at n points for
  \code{pcondgaussmixt} and density function for \code{dcondgaussmixt}.
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
\code{\link{condmixt.nll}}, \code{\link{condmixt.fwd}}
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

h <- 4 # number of hidden units
m <- 2 # number of components

# initialize a conditional mixture with Gaussian components and a dirac at zero
thetainit <- condgaussmixt.init(1,h,m,ytrain)

# compute mixture parameters 
params.mixt <- condgaussmixt.fwd(thetainit,h,m,t(xtrain))

cdf <- pcondgaussmixt(params.mixt,m,ytrain) # compute CDF 

}
