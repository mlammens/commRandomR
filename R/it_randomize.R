#' @export it_randomize
#' @title Randomize a matrix using the IT algorithm
#' @description
#' More info needed here
#'
#' @param matr A community matrix to be randomized.
#' @param n_iter Number of iterations to attempt to randomize the matrix.
#'
#' @return A randomized version of `matr` where the row
#' and column totals match that of `matr`, and the values
#' for each element were sampled proportional to the
#' corresponding element value in `matr`
#'
#' @details
#' Follows alogrithm as described in Ulrich and Gotelli 2010.
#' As per this source "IT assigns individuals randomly to matrix
#' cells with probabilities proportional to observed row and
#' column abundance totals until, for each row and column, total
#' abundances are reached."
#'
it_randomize <- function(matr, n_iter = 10000){

  # perform chisq test to get expected values
  #matr_exp <- chisq.test(matr)$expected
  matr_exp <- rowSums(matr) %o% colSums(matr)

  # make a new matrix of all zeros
  new_mat <- matr * 0

  # use a while loop to add matrices with a single, randomly
  # placed value of 1 until
  # the new_mat sum is the same as the original matrix
  iter <- 0
  while(sum(new_mat) < sum(matr)){
    # make a random matrix
    # this is done by randomly choosing a matrix index value, weighted
    # by the expected matrix values at each index
    zero_mat <- matr * 0
    zero_mat[sample(1:length(matr_exp), size = 1, prob = matr_exp)] <- 1

    # add the random matrix to the new_mat temporalily
    new_mat_temp <- new_mat + zero_mat

    # keep new_mat_temp if it doesn't result in over filling the
    # row or column sums
    if( !any((rowSums(new_mat_temp) - rowSums(matr)) > 0 ) &
        !any((colSums(new_mat_temp) - colSums(matr)) > 0 )) {
      new_mat <- new_mat_temp
    }

    # Re-adjust the weight matrix to remove rows or columns that are already filled
    rows2zero <- as.numeric(!(rowSums(new_mat) - rowSums(matr) == 0))
    matr_exp <- rows2zero * matr_exp

    # Increase the iter number
    iter <- iter + 1

    # Check if we've reached the max iterations
    if(iter >= n_iter){
      break()
    }
  }

  if(iter == n_iter){
    return(NA)
  } else {
    return(new_mat)
  }
}
