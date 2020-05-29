#' @export u_ratio_calculator
#' @title Calculate the U-ratio for a Community Matrix
#'
#' @description \code{u_ratio_calculator} returns the U-ratio
#' value for a community matrix, following the description of
#' this metric in Ulrich and Gotelli 2010 and Schroeder et al.
#' 2014, PLoS ONE.
#'
#' @param matr A community matrix with \emph{sites} as columns
#' and \emph{species} as rows
#'
#' @return U-ratio
#'
#' @details
#' As described in Schroeder et al. 2014, " U-ratio = V/W ; metric of co-occurrence, which compares the variance of row totals V with the sum of the column variances W. Low values of U indicate negative covariation in abundance between species. If species are segregated, we expect the U-ratio to be smaller than expected by chance.
#'
u_ratio_calculator <- function(matr){

  # Calculate variance of row totals (V)
  Var_Sum_Row <- var(rowSums(matr))

  # Calculate sum of column variances (W)
  Sum_Col_Var <- sum(apply(matr, MARGIN = 2, FUN = var))

  # U-Ratio = V/W
  U_Ratio <- (Var_Sum_Row / Sum_Col_Var)

  return(U_Ratio)
}
