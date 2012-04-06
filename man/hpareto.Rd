\name{hpareto}
\alias{hpareto}
\alias{dhpareto}
\alias{phpareto}
\alias{qhpareto}
\alias{rhpareto}
\title{
The Hybrid Pareto Distribution
}
\description{
Density, distribution function, quantile function and random generation
for the hybrid Pareto distribution with parameters xi, mu and sigma.
}
\usage{
dhpareto(y, xi, mu = 0, sigma = 1, log = FALSE, trunc = TRUE)
phpareto(q, xi, mu = 0, sigma = 1, trunc = TRUE)
qhpareto(p, xi, mu = 0, sigma = 1, trunc = TRUE)
rhpareto(n, xi, mu = 0, sigma = 1, trunc = TRUE)
}
\arguments{
  \item{y,q}{vector of quantiles}
  \item{p}{vector of probabilities}
  \item{n}{number of observations}
  \item{xi}{tail index parameter, inherited from the GPD}
  \item{mu}{location parameter, inherited from the Gaussian}
  \item{sigma}{scale parameter, inheristed from the Gaussian}
  \item{log}{logical, if TRUE, probabilities \code{p} are
    given as log(\code{p}).}
  \item{trunc}{logical, if TRUE (default), the hybrid Pareto density is
    truncated below zero.}
}
\details{
The hybrid Pareto density is given by a Gaussian with parameters
\code{mu} and \code{sigma} below the threshold
alpha (see the function \code{\link{hpareto.alpha}}) and by the GPD with
parameters \code{xi} and beta.(see the function
\code{\link{hpareto.beta}}) To
ensure continuity of the density and of its derivative at the threshold,
alpha and beta are appropriate functions of \code{xi}, \code{mu} and
\code{sigma}. Appropriate reweighting factor gamma ensures that the
density integrate to one.
}
\value{
\code{dhpareto} gives the density, \code{phpareto} gives the
distribution function, \code{qhpareto} gives the quantile function and
\code{rhpareto} generates random deviates.  
}
\references{
Carreau, J. and Bengio, Y. (2009), A Hybrid Pareto Model for
Asymmetric Fat-tailed Data: the Univariate Case, 12, Extremes
}
\author{
Julie Carreau
}

\seealso{
\code{\link{hpareto.alpha}}, \code{\link{hpareto.beta}} and
\code{\link{hpareto.gamma}}
}
