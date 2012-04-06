\name{lognormixt}
\alias{lognormixt}
\alias{dlognormixt}
\title{Mixture of Log-Normals}
\description{
Density  for a mixture of Log-Normals with
\code{m} components. 
}
\usage{
dlognormixt(params, x, log = FALSE)
}
\arguments{
  \item{params}{matrix of mixture parameters of dimension 3 x m, where m
    is the number of components, so that each column contains the
    mixture parameters (pi, mu, sigma) related to one component
}
  \item{x}{vector of sample data}
  \item{log}{logical, if TRUE, probabilities \code{p} are
    given as log(\code{p}).}
}

\value{\code{dlognormixt} gives the density of the mixture of Log-Normals}
\author{
Julie Carreau
}
\seealso{\code{\link{Lognormal}}}