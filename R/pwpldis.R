#' Discrete Piecewise Power-Law Distribution
#'
#' This set of functions provides the density, cumulative distribution function,
#' survival function, hazard function, quantile function, and random generation
#' for the discrete piecewise power-law distribution.
#'
#' @param x A numeric value to evaluate the function.
#' @param u A probability to evaluate the quantile function.
#' @param n An integer specifying the number of random samples to generate.
#' @param p A numeric vector specifying the minimum value \eqn{\tau_{(0)}} and
#' change points \eqn{\tau_{(1)}, \ldots, \tau_{(k)}} that define the partitions,
#' i.e., \eqn{p = (\tau_{(0)}, \tau_{(1)}, \ldots, \tau_{(k)})}. It must be sorted
#' in increasing order.
#' @param alpha A numeric vector of scaling parameters
#' \eqn{\alpha_{(1)}, \ldots, \alpha_{(k + 1)}} corresponding to each partition.
#' All values must be greater than 1.
#'
#' @return
#' - `dpwpldis()`: The probability mass function (PMF) evaluated at `x`.
#' - `ppwpldis()`: The cumulative distribution function (CDF) evaluated at `x`.
#' - `spwpldis()`: The survival function (1 - CDF) evaluated at `x`.
#' - `hpwpldis()`: The hazard function evaluated at `x`.
#' - `qpwpldis()`: The quantile function, returning the smallest `x` such that
#'                 \eqn{P(X \leq x) \geq u}
#' - `rpwpldis()`: A numeric vector of `n` random samples drawn from the distribution.
#'
#' @details The function assumes that `p` and `alpha` have the same length.
#'          Some functions rely on the Hurwitz zeta function from the
#'          VGAM package.
#'
#' @examples
#' p <- c(1, 3, 5)
#' alpha <- c(1.5, 2, 3)
#' dpwpldis(5, p, alpha)
#' ppwpldis(5, p, alpha)
#' spwpldis(5, p, alpha)
#' hpwpldis(5, p, alpha)
#' qpwpldis(0.5, p, alpha)
#' rpwpldis(10, p, alpha)
#'
#' @importFrom VGAM zeta
#' @importFrom stats runif
#'
#' @export
dpwpldis <- function (x, p, alpha)
{
  if (x < p[1]) stop(return(0))
  C <- Cj_pwpldis(p, alpha)
  j <- max(which(p <= x))

  x^(-alpha[j]) / zeta(alpha[j], shift = p[j]) * C[j]
}

#' @rdname dpwpldis
#' @export
ppwpldis <- function (x, p, alpha)
{
  if (x < p[1]) stop(return(0))
  C <- Cj_pwpldis(p, alpha)
  j <- max(which(p <= x))
  u <- 1 - zeta(alpha[j], shift = x + 1)/zeta(alpha[j], shift = p[j]) * C[j]

  return(u)
}

#' @rdname dpwpldis
#' @export
spwpldis <- function (x, p, alpha)
{
  1 - ppwpldis(x, p, alpha)
}

#' @rdname dpwpldis
#' @export
hpwpldis <- function (x, p, alpha)
{
  if (x < p[1]) stop(return(0))

  dpwpldis(x, p, alpha) / spwpldis(x - 1, p, alpha)
}

#' @rdname dpwpldis
#' @export
qpwpldis <- function (u, p, alpha)
{
  x <- p[1]
  while (TRUE) {
    if (ppwpldis(x, p, alpha) >= u) {
      break
    } else {
      x <- x + 1
    }
  }

  x
}

#' @rdname dpwpldis
#' @export
rpwpldis <- function(n, p, alpha)
{
  qq <- Vectorize(function(u) qpwpldis(u, p, alpha), "u")
  while(TRUE)
  {
    x <- qq(runif(n))
    nj <- n_each_interval(x, p)
    if (all(nj > 2)) break
  }

  x
}
