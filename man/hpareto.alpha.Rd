\name{hpareto.alpha}
\alias{hpareto.alpha}
\alias{hpareto.beta}
\alias{hpareto.gamma}
\title{
Auxillary Parameters of the Hybrid Pareto Distribution
}
\description{
Computes the junction point alpha, the GPD scale parameter beta and the
normalization factor gamma of the hybrid Pareto distributions.
}
\usage{
hpareto.alpha(xi, mu = 0, sigma = 1)
hpareto.beta(xi, sigma = 1)
hpareto.gamma(xi, mu = 0, sigma = 1, trunc = T)
}
\arguments{
  \item{xi}{Tail index parameter of the hybrid Pareto distribution.
}
  \item{mu}{Location parameter of the hybrid Pareto distribution.
}
  \item{sigma}{Scale parameter of the hybrid Pareto distribution.
}
\item{trunc}{Binary variable : if TRUE, the density of the hybrid Pareto
is truncated below zero}
}
\details{
Let z = (1+xi)^2/(2 pi) and W be the Lambert W function implemented in
\code{lambertw}. Then : 
alpha = mu + sigma * sqrt(W(z)), beta = sigma * |1 + xi| / sqrt(W(z))
and gamma is the integral from 0 (is trunc is true, -infinity otherwise)
to alpha. 
}
\value{
 Computation of the auxillary parameters alpha, beta and gamma.
}
\references{
Carreau, J. and Bengio, Y. (2009), A Hybrid Pareto Model for
Asymmetric Fat-tailed Data: the Univariate Case, 12, Extremes
}
\author{
Julie Carreau
}
\seealso{
\code{\link{dhpareto}},\code{\link{phpareto}} and \code{\link{rhpareto}}
}
\examples{
hpareto.alpha(0.1,0,1)
hpareto.beta(0.1,1)
hpareto.gamma(0.1,0,1)
}
