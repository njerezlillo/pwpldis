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
#' @param alpha A numeric vector of scaling parameters.
#' \eqn{\alpha_{(1)}, \ldots, \alpha_{(k + 1)}} corresponding to each partition.
#' All values must be greater than 1.
#'
#' @details
#' Consider a partition over \eqn{\mathcal{R} = [\tau_{(0)}, +\infty)} yielding
#' \eqn{\mathcal{R} = \mathcal{R}_1 \cup \ldots \cup \,\mathcal{R}_{k + 1}}, with
#' \eqn{\mathcal{R}_{j} = [\tau_{(j-1)}, \tau_{(j)})} for \eqn{j = 1, \dots, k + 1},
#' where \eqn{\tau_{(0)} \geq 1} is the sample minimum, \eqn{\tau_{(1)} < \ldots < \tau_{(k)}}
#' are change points, and \eqn{\tau_{(k+1)} = +\infty}. The probability mass function of
#' the discrete piecewise power-law model is given by:
#' \deqn{p(x)=\displaystyle\sum_{j=1}^{k + 1} \left[\dfrac{x^{-\alpha_j}}{\zeta(\alpha_j,\tau_{(j-1)})}
#' \cdot C_{j-1} \right] \text{I}_{\mathcal{R}_j}(x)}
#' where \eqn{C_0=1}, and \eqn{C_{i}=\displaystyle\prod_{h=1}^{i} \frac{\zeta(\alpha_h,\tau_{(h)})}{\zeta(\alpha_{h},\tau_{(h-1)})},}
#' with \eqn{i=1,\ldots,k}, normalization constants.
#'
#' The function assumes that `p` and `alpha` have the same length.
#'
#' Some functions rely on the Hurwitz zeta function from the VGAM package.
#'
#' This function was inspired by the implementation of the `PWEXP`
#' package in `R`. Visit their GitHub repository
#' [PWEXP](https://github.com/zjph602xtc/PWEXP).
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
#' @references
#' Jerez-Lillo, N., Rodrigues, F. A., Ferreira, P. H., & Ramos, P. L. (2025).
#' Beyond the Power Law: Estimation, Goodness-of-Fit, and a Semiparametric
#' Extension in Complex Networks. arXiv preprint arXiv:2311.11200. Available at:
#' \url{https://arxiv.org/abs/2311.11200}
#'
#' @examples
#' # Define parameters for the piecewise discrete power-law distribution
#' p <- c(1, 3, 5)
#' alpha <- c(1.5, 2, 3)
#'
#' # Calculate the density at a specific point
#' dpwpldis(5, p, alpha)
#'
#' # Calculate the cumulative distribution function at point 5
#' ppwpldis(5, p, alpha)
#'
#' # Calculate the survival function at point 5
#' spwpldis(5, p, alpha)
#'
#' # Calculate the hazard function at point 5
#' hpwpldis(5, p, alpha)
#'
#' # Calculate the quantile corresponding to a 50% probability (median)
#' qpwpldis(0.5, p, alpha)
#'
#' # Generate 10 random samples from the discrete power-law distribution
#' rpwpldis(10, p, alpha)
#'
#' @seealso [Cj_pwpldis]
#'
#' @importFrom VGAM zeta
#' @importFrom stats runif
#'
#' @export
dpwpldis <- Vectorize(function(x, p, alpha) {
  if (x < p[1]) return(0)
  C <- Cj_pwpldis(p, alpha)
  j <- max(which(p <= x))

  x^(-alpha[j]) / zeta(alpha[j], shift = p[j]) * C[j]
}, vectorize.args = "x")

#' @rdname dpwpldis
#' @export
ppwpldis <- Vectorize(function(x, p, alpha) {
  if (x < p[1]) return(0)
  C <- Cj_pwpldis(p, alpha)
  j <- max(which(p <= x))

  1 - zeta(alpha[j], shift = x + 1) / zeta(alpha[j], shift = p[j]) * C[j]
}, vectorize.args = "x")

#' @rdname dpwpldis
#' @export
spwpldis <- Vectorize(function(x, p, alpha) {
  1 - ppwpldis(x, p, alpha)
}, vectorize.args = "x")

#' @rdname dpwpldis
#' @export
hpwpldis <- Vectorize(function(x, p, alpha) {
  if (x < p[1]) return(0)
  dpwpldis(x, p, alpha) / spwpldis(x - 1, p, alpha)
}, vectorize.args = "x")

#' @rdname dpwpldis
#' @export
qpwpldis <- Vectorize(function(u, p, alpha) {
  x <- p[1]
  while (TRUE) {
    if (ppwpldis(x, p, alpha) >= u) {
      break
    } else {
      x <- x + 1
    }
  }
  x
}, vectorize.args = "u")

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
