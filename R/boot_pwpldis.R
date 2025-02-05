#' Bootstrap for Piecewise Power-Law Distribution
#'
#' This function performs a bootstrap procedure to estimate the distribution of
#' parameters for a piecewise power-law model. It resamples the dataset multiple
#' times and fits the model parameters for each bootstrap sample.
#'
#' @param time A numeric vector of observations.
#' @param nsim An integer specifying the number of bootstrap simulations.
#' Default is 100.
#' @param brks A numeric vector of change points for the piecewise power-law model.
#' If NULL, the change points are estimated.
#' @param nbreak An integer specifying the number of change points in the model.
#' Default is 1.
#' @param exclude_int A numeric vector of two values specifying an interval where
#' change points should not be placed (e.g., `exclude_int = c(5, Inf)` excludes
#' change points greater than 5).
#' @param min_pt_tail An integer specifying the minimum number of observations
#' required in the last segment to ensure stable parameter estimation. Default is 5.
#' @param max_set An integer specifying the maximum number of possible change point
#' combinations considered for model fitting. Default is 1000.
#' @param tol A numeric value specifying the tolerance for combining close change points.
#' Default is 1e-4.
#' @param parallel A logical value indicating whether to use parallel processing
#' for bootstrap simulations. Default is FALSE.
#' @param mc.core An integer specifying the number of cores to use in parallel processing.
#' Default is 4.
#' @param ... Additional arguments passed to the underlying `fit_pwpldis` function.
#'
#' @details
#' This function performs a bootstrap procedure by resampling the original dataset
#' (`time`) with replacement for `nsim` iterations. For each resample, the function
#' fits a piecewise power-law model using either the specified or estimated change points,
#' and the model parameters are recorded.
#'
#' The `exclude_int` parameter prevents change points from being placed in a specific interval.
#' The `min_pt_tail` parameter ensures that the last segment has a sufficient number of
#' observations for stable parameter estimation.
#'
#' After the bootstrap simulations, the results are aggregated and returned,
#' including the estimated change points and scaling parameters for each segment.
#'
#' This function was inspired by the implementation of the `PWEXP`
#' package in `R`. Visit their GitHub repository
#' [PWEXP](https://github.com/zjph602xtc/PWEXP).
#'
#' @return A `data.frame` with the results of the bootstrap procedure, including the
#' following columns:
#' \describe{
#'   \item{taus}{Estimated change points from each bootstrap iteration.}
#'   \item{alphas}{Estimated scaling parameters for each segment from each
#'   bootstrap iteration.}
#'   \item{likelihood}{The log-likelihood in each bootstrap iteration.}
#'   \item{AIC}{The Akaike Information Criterion in each bootstrap iteration.}
#'   \item{BIC}{The Bayesian Information Criterion in each bootstrap iteration.}
#' }
#'
#' @references
#' Jerez-Lillo, N., Rodrigues, F. A., Ferreira, P. H., & Ramos, P. L. (2025).
#' Beyond the Power Law: Estimation, Goodness-of-Fit, and a Semiparametric
#' Extension in Complex Networks. arXiv preprint arXiv:2311.11200. Available at:
#' \url{https://arxiv.org/abs/2311.11200}
#'
#' @examples
#'
#' @seealso [fit_pwpldis]
#'
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach `%dopar%`
#' @importFrom parallel makeCluster stopCluster
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
boot_pwpldis <- function(time, nsim = 100, brks = NULL,
                         nbreak = 1, exclude_int = NULL,
                         min_pt_tail = 5, max_set = 1000,
                         tol = 1e-4, parallel = FALSE,
                         mc.core = 4, ...) {

  dat <- data.frame(time=time)
  n <- NROW(dat)
  res_all <- fit_pwpldis(time=dat$time, breakpoint=brks, nbreak=nbreak, exclude_int=exclude_int, min_pt_tail=min_pt_tail, max_set=max_set, trace=FALSE, tol=tol)

  ind <- order(dat$time)
  dat <- dat[ind,,drop=F]
  if (tol!=0){
    dat$time <- combine_time(dat$time, tol)
  }

  if (nbreak==0){
    nbreak <- length(attr(res_all, 'lam'))-1
  }
  pb <- txtProgressBar(max = nsim, style = 3) # optional

  if (parallel){
    doSNOW::registerDoSNOW(cl <- parallel::makeCluster(mc.core))
    `%dopar%` <- foreach::`%dopar%`
    res_all_tp <- foreach::foreach(i=1:(nsim-1), .combine = 'rbind', .inorder = FALSE, .errorhandling = 'remove', .packages = 'PWEXP', .options.snow=list(progress=function(n)setTxtProgressBar(pb, n))) %dopar% {
      dat_b <- dat[sample.int(n, n, replace = T), 1]
      res <- suppressWarnings(fit_pwpldis(time=dat_b, breakpoint=brks, nbreak=nbreak, exclude_int=exclude_int, min_pt_tail=min_pt_tail, max_set=max_set, trace=FALSE, tol=0))
    }
    res_all <- rbind(res_all, res_all_tp)
    parallel::stopCluster(cl)
  }else{
    for (i in 1:(nsim-1)){
      setTxtProgressBar(pb, i) # optional
      dat_b <- dat[sample.int(n, n, replace = T), 1]
      if (all(n_each_interval(dat_b, c(min(time), brks)) > 1)) {
        res <- suppressWarnings(fit_pwpldis(time=dat_b, breakpoint=brks, nbreak=nbreak, exclude_int=exclude_int, min_pt_tail=min_pt_tail, max_set=max_set, trace=FALSE, tol=0))
        res_all <- rbind(res_all, res)
      }
    }
  }

  setTxtProgressBar(pb, nsim) # optional
  close(pb) # optional
  res_all[is.infinite(res_all[,1]),] <- NA
  res_all[is.na(res_all[,1]),] <- suppressWarnings(matrix(colMeans(res_all, na.rm=T), ncol=NCOL(res_all), nrow=sum(is.na(res_all[,1])), byrow = T))
  if (nbreak!=0){
    attr(res_all,'brk') <- res_all[,1:nbreak,drop=F]
  }else{
    attr(res_all,'brk') <- NULL
  }
  attr(res_all,'lam') <- res_all[,(nbreak+1):(2*nbreak+1),drop=F]
  class(res_all) <- c('fit.boot.pwpldis', 'data.frame')
  return(res_all)
}
