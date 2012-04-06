\name{softplus}
\alias{softplus}
\alias{softplusinv}
\title{
  Softplus Transform
}
\description{
The softplus (and inverse softplus) transform is useful to introduce
positivity constraints on parameters of a function that will be
optimized (e.g. MLE of the scale parameter of a density function). The
softplus has been introduced to replace the exponential which might blow
up for large argument. The softplus is given by : log(1+exp(x)) and
converges to x for large values of x. Some care has been taken in the
implementation of the softplus function to handle some numerical
issues.
}
\usage{
softplus(x)
softplusinv(y)
}
\arguments{
  \item{x}{ is the value of the unconstrained parameter which is optimized    
  }
  \item{y}{ is the value of the positively constrained parameter}}
\details{
Let sigma be the scale parameter of a density for which maximumm
likelihood estimation will be performed. Then we can consider optimizing
\code{softplusinv}(sigma) to ensure positivity of this parameter. Let sigma.unc
be the optimzed unconstrained parameter, then \code{softplus}(sigma.unc) is the
value of the MLE.
}
\value{
The value of the softplus (or inverse softplus) transform.
}
\references{
Dugas, C., Bengio, Y., Belisle, F., Nadeau, C. and Garcia,
R. (2001), A universal approximator of convex functions applied to
option pricing, 13, Advances in Neural Information Processing Systems
}
\author{
Julie Carreau
}

