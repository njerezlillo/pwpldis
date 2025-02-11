#' Combine Close Observations Based on a Tolerance Threshold
#'
#' This function modifies a numeric vector of observations by merging values
#' that are within a specified tolerance threshold.
#'
#' @param time A numeric vector containing the observations to be adjusted.
#' @param tol A numeric value representing the tolerance threshold.
#'
#' @return A modified numeric vector where values within the tolerance
#' threshold have been merged with the preceding value.
#'
#' @author Tianchen Xu (<zjph602xutianchen@gmail.com>)
#'
#' @references This function was extracted from the GitHub repository
#' [PWEXP](https://github.com/zjph602xtc/PWEXP).
#'
#' @export
combine_time <- function(time, tol){
  f_tmp <- time[1]
  f_N <- length(time)
  tol <- (time[f_N] - f_tmp) * tol
  for (i in 2:f_N){
    val <- time[i]
    if ((val - f_tmp) < tol){
      time[i] <- f_tmp
    }else{
      f_tmp <- val
    }
  }
  return(time)
}
