#' Log-Likelihood of the Discrete Piecewise Power-Law Model
#'
#' This function computes the log-likelihood for a discrete piecewise power-law model.
#' It evaluates the log-likelihood given a set of parameters, observed data, and
#' the change points that define the intervals.
#'
#' @param alpha A numeric vector of scaling parameters
#' \eqn{\alpha_{(1)}, \ldots, \alpha_{(k + 1)}} corresponding to each partition.
#' All values must be greater than 1.
#' @param x A numeric vector of observed data values.
#' @param p A numeric vector specifying the minimum value \eqn{\tau_{(0)}} and
#' change points \eqn{\tau_{(1)}, \ldots, \tau_{(k)}} that define the partitions,
#' i.e., \eqn{p = (\tau_{(0)}, \tau_{(1)}, \ldots, \tau_{(k)})}. It must be sorted
#' in increasing order.
#'
#' @return A numeric value representing the log-likelihood of the data given
#' the specified parameters and change points
#'
#' @details The function assumes that `p` and `alpha` have the same length.
#'          The calculation involves the Hurwitz zeta function from the
#'          VGAM package.
#'
#' @examples
#' x <- c(7, 7, 6, 4, 1, 5, 2, 6, 3, 5)
#' loglik_pwpldis(c(2.0, 3.5), x, c(1, 4))
#' loglik_pwpldis(c(2.5, 3.0), x, c(1, 4))
#' loglik_pwpldis(c(1.5, 2.0), x, c(1, 4))
#'
#' @importFrom VGAM zeta
#'
#' @export
loglik_pwpldis <- function(alpha, x, p)
{
  k <- length(p)
  index <- index_each_interval(x, p)
  aux <- vector(length = k)

  for (j in 1:k) {
    aux[j] <- sum(alpha[j] * log(x[index[[j]]]) +
                    log(zeta(alpha[j], shift = p[j])))
  }

  logC <- log(Cj_pwpldis(p, alpha))[-(length(p) + 1)]
  nj <- n_each_interval(x, p)

  sum(nj * logC - aux)
}
