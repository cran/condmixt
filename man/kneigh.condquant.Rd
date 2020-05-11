\name{kneigh.condquant}
\alias{kneigh.condquant}
\title{
  Conditional quantile estimation from nearest neighbors.
  }
\description{
Conditional quantile estimation from k-nearest neighbors in the
explanatory variable space.
}
\usage{
kneigh.condquant(x, y, k = 10, p = 0.9)
}
\arguments{
  \item{x}{
    Matrix of explanatory (independent) variables of dimension d x n, d is
    the number of variables and n is the number of examples (patterns)
}
\item{y}{
  Vector of n dependent variables
}
  \item{k}{
    Number of neighbors, default is 10.
}
  \item{p}{
    Probability level, default is 0.99.
}
}
\details{
For each example j (each column) in the matrix \code{x}, its \code{k}
nearest neighbors in terms of Euclidean distance are identified. Let
j1,..., jk be the \code{k} nearest neighbors. Then, the conditional
quantile is estimated by computing the sample quantile over
y[j1],...,y[jk].
}
\value{
A vector of quantile of length n.  
}
\references{
Bishop, C. (1995), Neural Networks for Pattern Recognition, Oxford
}
\author{
Julie Carreau
}

\seealso{
 \code{\link{quantile}}
}
\examples{
# generate train data
ntrain <- 500
xtrain <- runif(ntrain) 
ytrain <- rfrechet(ntrain,loc = 3*xtrain+1,scale =
0.5*xtrain+0.001,shape=xtrain+1)
plot(xtrain,ytrain,pch=22) # plot train data
qgen <- qfrechet(0.99,loc = 3*xtrain+1,scale =
 0.5*xtrain+0.001,shape=xtrain+1) # compute quantile from generative model
points(xtrain,qgen,pch=".",col="orange")

kquant <- kneigh.condquant(t(xtrain),ytrain,p=0.99) # compute estimated quantile
                                                               
points(xtrain,kquant,pch="o",col="blue")
# sample quantiles are not good in the presence of heavy-tailed data


ytrain <- rlnorm(ntrain,meanlog = 3*xtrain+1,sdlog =
0.5*xtrain+0.001)
dev.new()
plot(xtrain,ytrain,pch=22) # plot train data
qgen <- qlnorm(0.99,meanlog = 3*xtrain+1,sdlog =
 0.5*xtrain+0.001) # compute quantile from generative model
points(xtrain,qgen,pch=".",col="orange")
# compute estimated quantile
kquant <- kneigh.condquant(t(xtrain),ytrain,p=0.99) 

points(xtrain,kquant,pch="o",col="blue") # a bit more reasonable 

}
