\name{gaussmixt}
\alias{gaussmixt}
\alias{dgaussmixt}
\alias{pgaussmixt}
\title{Mixture of Gaussians}
\description{
Density and distribution function for a mixture of Gaussians with
\code{m} components. 
}
\usage{
dgaussmixt(params, x, log = FALSE, trunc = TRUE)
pgaussmixt(params, x, trunc = TRUE)
}
\arguments{
  \item{params}{matrix of mixture parameters of dimension 3 x m, where m
    is the number of components, so that each column contains the
    mixture parameters (pi, mu, sigma) related to one component
}
  \item{x}{vector of sample data}
  \item{log}{logical, if TRUE, probabilities \code{p} are
    given as log(\code{p}).}
  \item{trunc}{logical, if TRUE (default), the Gaussian density is
    truncated below zero.}
}

\value{\code{dgaussmixt} gives the density and \code{pgaussmixt} gives the
distribution function}
\author{
Julie Carreau
}
\seealso{\code{\link{Normal}}}