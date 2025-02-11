#' Fit the Discrete Piecewise Power-law Distribution
#'
#' This function fits a piecewise power-law distribution to datasets. The user
#' can specify all change points or allow the function to estimate them.
#' It also provides the option to exclude intervals and specify
#' the minimum number of observations required to estimate the last parameter.
#'
#' @param time A numeric vector of observations.
#' @param breakpoint A numeric vector specifying fixed change points.
#' @param nbreak An integer specifying the number of change points in the model.
#' @param exclude_int A numeric vector of two values defining an interval in which
#' estimated change points should be excluded (e.g., `exclude_int = c(5, Inf)`
#' excludes change points after `time = 5`).
#' @param min_pt_tail An integer specifying the minimum number of observations
#' required to estimate the last parameter.
#' @param max_set An integer specifying the maximum number of possible change points
#' combinations to consider.
#' @param trace A logical value indicating whether to display the optimization
#' process. Default is `FALSE`.
#' @param tol A numeric value specifying the minimum allowed gap between two
#' change points.
#'
#' @details
#' If the user specifies change points, the function checks their validity and removes
#' change points after the maximum or before the minimum. The `exclude_int` argument
#' allows excluding change points that are too close to the tail. The `min_pt_tail`
#' argument helps in estimating a more robust scaling parameter for the tail by
#' requiring a minimum number of observations.
#'
#' This function was inspired by the implementation of the `PWEXP`
#' package in `R`. Visit their GitHub repository
#' [PWEXP](https://github.com/zjph602xtc/PWEXP).
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{taus}{Estimated change points.}
#'   \item{alphas}{Estimated scaling parameters for each piece.}
#'   \item{likelihood}{The log-likelihood of the model.}
#'   \item{AIC}{The Akaike Information Criterion for the model.}
#'   \item{BIC}{The Bayesian Information Criterion for the model.}
#' }
#'
#' @references
#' Jerez-Lillo, N., Rodrigues, F. A., Ferreira, P. H., & Ramos, P. L. (2025).
#' Beyond the Power Law: Estimation, Goodness-of-Fit, and a Semiparametric
#' Extension in Complex Networks. arXiv preprint arXiv:2311.11200. Available at:
#' \url{https://arxiv.org/abs/2311.11200}
#'
#' @examples
#' # Initial specification
#' set.seed(2025)
#' p <- c(1, 3, 5)
#' alpha <- c(1.5, 2, 3)
#'
#' # Generate a dataset from the discrete piecewise power-law model
#' df <- rpwpldis(300, p, alpha)
#'
#' # Fit the model with fixed change points (cps) at 3 and 5
#' fit_pwpldis(df, breakpoint = c(3, 5))
#'
#' # Fit the model allowing it to estimate 1 cp
#' fit_pwpldis(df, nbreak = 1)
#'
#' # Fit the model with 1 estimated cp and display optimization process
#' fit_pwpldis(df, nbreak = 1, trace = TRUE)
#'
#' # Fit the model with 1 estimated cp, requiring at least 100 observations in the last partition
#' fit_pwpldis(df, nbreak = 1, min_pt_tail = 100)
#'
#' # Compute the number of observations in each interval defined by cps
#' n_each_interval(df, c(1, 4))
#'
#' # Fit the model allowing it to estimate 2 cps
#' fit_pwpldis(df, nbreak = 2)
#'
#' # Fit the model with 2 estimated cps, excluding cps beyond 5
#' fit_pwpldis(df, nbreak = 2, exclude_int = c(5, Inf))
#'
#' @seealso [dpwpldis]
#'
#' @importFrom fastmatch ctapply
#' @importFrom maxLik maxLik
#'
#' @export
fit_pwpldis <- function(time, breakpoint = NULL, nbreak = NULL,
                        exclude_int = NULL, min_pt_tail = 2,
                        max_set = 10000, trace = FALSE, tol = 1e-4){

  time <- sort(time)
  breakpoint <- sort(breakpoint)

  N <- length(time)
  tau_0 <- min(time)
  optimizer <- "mle"
  n_fix_brk <- length(breakpoint)

  if (!is.null(breakpoint)) nbreak <- n_fix_brk

  if (n_fix_brk==0  && is.null(nbreak)){
    nbreak <- ceiling(8 * length(time)^0.2)
    message(paste0('Number of change points = ', nbreak))
  }

  while (n_fix_brk > 0){
    tmpi <- findInterval(time, vec=c(tau_0, breakpoint, Inf))
    numerator <- ctapply(rep(1, N), tmpi, sum) # ctapply(event, tmpi, sum) #EVENTO YA NO EXISTE
    if (all(numerator > 1) && all((1:(n_fix_brk+1)) %in% names(numerator))){
      break
    }
    ind <- min(which(numerator <= 1), setdiff(1:(n_fix_brk+1), names(numerator)), na.rm=T)
    if (is.na(breakpoint[ind])){
      warning(paste0(breakpoint[ind-1],' is too large. No event after this breakpoint. We will remove it.'))
      breakpoint <- breakpoint[-(ind-1)]
    }else if (ind==1){
      warning(paste0(breakpoint[ind],' is too early. No event before this breakpoint. We will remove it.'))
      breakpoint <- breakpoint[-ind]
    }else{
      warning(paste0('There are NO observations between ', breakpoint[ind-1], ' and ',
                     breakpoint[ind], '. We will use their mid point ',
                     mean(breakpoint[(ind-1):ind]), ' as a breakpoint. '))
      breakpoint[ind-1] <- mean(breakpoint[(ind-1):ind])
      breakpoint <- breakpoint[-ind]
    }
    n_fix_brk <- length(breakpoint)
    if (n_fix_brk == 0){
      breakpoint <- NULL
    }
  }

  if (!is.null(breakpoint)){
    if (nbreak==0 | nbreak == length(breakpoint)){
      breakpoint <- matrix(breakpoint, nrow=1)
      optimizer <- "fixed"
    }
    if (nbreak!=0 & nbreak < length(breakpoint)){
      stop('nbreak is the total number of change points, which must be equal or larger than the length of \'breakpoint\'')
    }
  }

  if (optimizer == 'mle') {
    change_cand <- FALSE
    if(min_pt_tail!=0 | min_pt_tail!=1){
      change_cand <- TRUE
      if ((length(time)-min_pt_tail) < 2){
        warning('Insufficient data due to too many points reserved for tail. Ignore \'min_pt_tail\'.')
        candidate_time <- time
      }else{
        candidate_time <- time[time <= time[length(time)-min_pt_tail]]
        if (length(candidate_time) <= 2*(nbreak-length(breakpoint))){
          warning('Insufficient data due to too many points reserved for tail. Ignore \'min_pt_tail\'.')
          candidate_time <- time
        }
      }
    }
    if(!is.null(exclude_int)){
      change_cand <- TRUE
      candidate_time <- time[time < exclude_int[1] | time > exclude_int[2]]
      if (length(candidate_time) <= 2*(nbreak-length(breakpoint))){
        warning('Insufficient data due to wide coverage of \'exclude_int\'. Ignore \'exclude_int\'.')
        candidate_time <- time
      }
    }
    if (!change_cand){
      candidate_time <- time
    }
    breakpoint <- get_grid(candidate_time, breakpoint, nbreak, max_set)
  }

  if (is.infinite(breakpoint[1, 1])){
    res <- matrix(-Inf, ncol = 2 * NCOL(breakpoint) + 3, nrow = 1)
  } else {
    res <- matrix(-Inf, ncol = 2 * NCOL(breakpoint) + 3, nrow = NROW(breakpoint))
    for (i in 1:NROW(breakpoint)){
      brk0 <- breakpoint[i,]
      brk <- c(tau_0, brk0)
      numerator <- n_each_interval(time, brk)
      if (any(numerator <= 1)){
        next
      }
      l <- function(z) loglik_pwpldis(z, time, brk)
      f <- maxLik(l, start = runif(length(brk), 1.5, 3.5),
                  constraints = list(ineqA = diag(rep(1, length(brk))),
                                     ineqB = rep(-1.01, length(brk))))

      loglikelihood <- f$maximum
      res[i,] <- c(brk, f$estimate, loglikelihood)
    }
  }

  res <- data.frame(res)
  if (!trace){
    res <- res[which.max(res[,NCOL(res)])[1],,drop=F]
  }
  n_k <- 2 * max(n_fix_brk, nbreak) + 1
  aic <- 2 * (n_k-n_fix_brk) - 2 * res[,NCOL(res)]
  bic <- (n_k - n_fix_brk) * log(N) - 2 * res[,NCOL(res)]
  res <- cbind(res, aic, bic)
  if (!trace){
    attr(res,'alpha') <- as.numeric(res[,(NCOL(breakpoint)+1):(2*NCOL(breakpoint)+1)])
    attr(res,'tau_') <- as.numeric(res[,1:NCOL(breakpoint)])
  }
  colnames(res) <- c(paste0('tau_', 0:NCOL(breakpoint)), paste0('alpha', 1:(NCOL(breakpoint)+1)), 'likelihood','AIC','BIC')
  attr(res,'para') <- list(time = time, breakpoint = breakpoint, nbreak = nbreak, exclude_int=exclude_int, min_pt_tail=min_pt_tail)

  return(res)
}
