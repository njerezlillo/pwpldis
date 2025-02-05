#' Generate a Grid of Breakpoints
#'
#' This function constructs a grid of change points for the discrete piecewise
#' power-law model. It selects a subset of unique observations and ensures that
#' the number of combinations does not exceed a specified limit.
#'
#' @param time A numeric vector of observations.
#' @param breakpoint A numeric vector specifying fixed change points.
#' @param nbreak An integer specifying the number of change points in the model.
#' @param max_set An integer specifying the maximum number of combinations to generate.
#' Default is 5000.
#' @param remove_first A logical value indicating whether to exclude the first time point.
#' Default is `TRUE`.
#'
#' @return A matrix where each row represents a combination of observations
#' forming a set of change points for the model.
#'
#' @author Tianchen Xu (<zjph602xutianchen@gmail.com>)
#'
#' @references This function was extracted from the GitHub repository
#' [PWEXP](https://github.com/zjph602xtc/PWEXP).
#'
#' @importFrom utils combn
#'
#' @export
get_grid <- function(time, breakpoint, nbreak, max_set = 5000, remove_first = TRUE){
  time <- unique(time)
  if (remove_first){
    time <- time[-1]
  }

  if (nbreak==0){
    nbreak <- length(breakpoint)
  }
  if (!is.null(breakpoint)){
    nbreak <- nbreak - length(breakpoint)
    time <- setdiff(time, breakpoint)
    except <- sapply(breakpoint, function(x) findInterval(x, time))
    time <- time[-c(except, except+1)]
  }
  nr <- length(time)
  nl <- 1

  if (choose(nr, nbreak) > max_set){
    while((nr-nl)>1.1){
      if (choose(floor((nl+nr)/2), nbreak) > max_set){
        nr <- floor((nl+nr)/2)
      }else{
        nl <- floor((nl+nr)/2)
      }
    }
  }
  time <- sort(sample(time,nr,replace = F))
  if (length(time) < nbreak){
    setting <- matrix(-Inf, ncol=nbreak, nrow=1)
  }else{
    setting <- t(combn(time, nbreak))
  }


  if (NROW(setting) > max_set){
    setting <- setting[sort(sample.int(NROW(setting), max_set, replace = F)),,drop=F]
  }

  if (!is.null(breakpoint)){
    setting <- t(apply(setting, 1, function(x)sort(c(breakpoint, x))))
  }
  return(setting)
}
