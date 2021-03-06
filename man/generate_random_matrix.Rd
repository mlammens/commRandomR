% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_random_matrix.R
\name{generate_random_matrix}
\alias{generate_random_matrix}
\title{Generate a Random Matrix}
\usage{
generate_random_matrix()
}
\value{
A new random matrix representing an ecological community with
sites and species and associated abundances.
}
\description{
\code{generate_random_matrix} returns a community matrix
generated following the guidelines in Ulrich and Gotelli 2010.
}
\details{
This function generates a random community matrix following the
algorithm documented in Ulrich and Gotelli 2010 for the \eqn{M_R}
type matrix.

In the resulting matrix, each column represents a different \emph{site}
and each row represents a different \emph{species}. The matrix
elements represent the \emph{abundance} of a specific species
at a specific site.
}
