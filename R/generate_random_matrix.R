#' @export generate_random_matrix
#' @title Generate a Random Matrix
#'
#' @description \code{generate_random_matrix} returns a community matrix
#' generated following the guidelines in Ulrich and Gotelli 2010.
#'
#' @return A new random matrix representing an ecological community with
#' sites and species and associated abundances.
#'
#' @details
#' This function generates a random community matrix following the
#' algorithm documented in Ulrich and Gotelli 2010 for the \eqn{M_R}
#' type matrix.
#'
#' In the resulting matrix, each column represents a different \emph{site}
#' and each row represents a different \emph{species}. The matrix
#' elements represent the \emph{abundance} of a specific species
#' at a specific site.
#'
generate_random_matrix <- function(){

  # Choose number of columns (sites)
  n_sites <- sample(5:50, size = 1)

  # Choose number of species in the meta-community
  n_spec_metacomm <- sample(10:200, size = 1)

  # Determine the total abundance N_i of species i
  # using N_i = exp(x_i/2a), where x_i is ~ N(0,1)
  # and a is Unif(0.1,1)
  spec_abund <- exp(rnorm(n_spec_metacomm) / (2*runif(1, min = 0.1, max = 1)) )

  # Sort the spec_abund values
  spec_abund_sort <- sort(spec_abund, decreasing = TRUE)

  # Generate a "max" cutoff point
  max_cutoff <- sample((n_spec_metacomm/2):n_spec_metacomm, size = 1)

  # Get rleative species abundances
  spec_relabund <- spec_abund_sort[1:max_cutoff]

  # Generate relative site carrying capacities
  site_k <- runif(n = n_sites, min = 0, max = 1)

  # Calculate matrix where each element is the row column product
  wt_matrix <- spec_relabund %o% site_k

  ## Create a new matrix using the random fill approach

  # make a new matrix of all zeros
  new_mat <- wt_matrix * 0

  # use a while loop to add matrices with a single, randomly
  # placed value of 1 until
  # the new_mat sum is the same as the original matrix
  while(any(rowSums(new_mat) == 0) | any(colSums(new_mat) == 0)){
    # make a random matrix
    # this is done by randomly choosing a matrix index value, weighted
    # by the expected matrix values at each index
    zero_mat <- wt_matrix * 0
    zero_mat[sample(1:length(wt_matrix), size = 1, prob = wt_matrix)] <- 1

    # add the random matrix to the new_mat
    new_mat <- new_mat + zero_mat
  }


  return(new_mat)
}
