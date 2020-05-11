\name{hparetomixt}
\alias{hparetomixt}
\alias{dhparetomixt}
\alias{phparetomixt}
\title{Mixture of hybrid Paretos}
\description{
Density and distribution function for a mixture of hybrid Paretos with
\code{m} components. 
}
\usage{
dhparetomixt(params, x, log = FALSE, trunc = TRUE)
phparetomixt(params, x, trunc = TRUE)
}
\arguments{
  \item{params}{matrix of mixture parameters of dimension 4 x m, where m
    is the number of components, so that each column contains the
    mixture parameters (pi, xi, mu, sigma) related to one component
}
  \item{x}{vector of sample data}
  \item{log}{logical, if TRUE, probabilities \code{p} are
    given as log(\code{p}).}
  \item{trunc}{logical, if TRUE (default), the hybrid Pareto density is
    truncated below zero.}
}

\value{\code{dhparetomixt} gives the density and \code{phparetomixt} gives the
distribution function}
\references{
Carreau, J. and Bengio, Y. (2009), A Hybrid Pareto Model for
Asymmetric Fat-tailed Data: the Univariate Case, 12, Extremes
}
\author{
Julie Carreau
}
\seealso{\code{\link{hparetomixt.init}}}