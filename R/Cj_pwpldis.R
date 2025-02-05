#' Normalizing Constants for the Discrete Piecewise Power-Law Model
#'
#' This function computes the normalizing constants required for the discrete
#' piecewise power-law model, ensuring that the probability mass function sums
#' to one over the defined partitions.
#'
#' @param p A numeric vector specifying the minimum value \eqn{\tau_{(0)}} and
#' change points \eqn{\tau_{(1)}, \ldots, \tau_{(k)}} that define the partitions,
#' i.e., \eqn{p = (\tau_{(0)}, \tau_{(1)}, \ldots, \tau_{(k)})}. It must be sorted
#' in increasing order.
#' @param alpha A numeric vector of scaling parameters
#' \eqn{\alpha_{(1)}, \ldots, \alpha_{(k + 1)}} corresponding to each partition.
#' All values must be greater than 1.
#'
#' @details
#' The constants are given by:
#' \deqn{C_0=1, \quad C_{i}=\displaystyle\prod_{h=1}^{i} \frac{\zeta(\alpha_h,\tau_{(h)})}{\zeta(\alpha_{h},\tau_{(h-1)})},}
#' with \eqn{i=1,\ldots,k.}
#'
#' The function assumes that `p` and `alpha` have the same length.
#'
#' The calculation involves the Hurwitz zeta function from the VGAM package.
#'
#' @return A numeric vector containing the normalizing constants for each
#'         partition of the model. By definition, the first value returned
#'         is 1, and the last is 0.
#'
#' @references
#' Jerez-Lillo, N., Rodrigues, F. A., Ferreira, P. H., & Ramos, P. L. (2025).
#' Beyond the Power Law: Estimation, Goodness-of-Fit, and a Semiparametric
#' Extension in Complex Networks. arXiv preprint arXiv:2311.11200. Available at:
#' \url{https://arxiv.org/abs/2311.11200}
#'
#' @examples
#' # Example with increasing parameters
#' Cj_pwpldis(c(1, 3, 5), c(1.5, 2.0, 2.5))
#'
#' # Example with decreasing parameters
#' Cj_pwpldis(c(1, 3, 5), c(2.5, 2.0, 1.5))
#'
#' # Example with non-equidistant partitions
#' Cj_pwpldis(c(1, 2, 6), c(2.0, 3.0, 1.5))
#'
#' @importFrom VGAM zeta
#'
#' @export
Cj_pwpldis <- function (p, alpha)
{
  if (any(sort(p) != p)) {
    stop(paste("p must be crecient", "\n", ""))
  }
  if (any(alpha <= 1)) {
    stop(paste("alpha must be > 1", "\n", ""))
  }

  aux <- c(1, rep(NA, length(p) - 1), 0)

  for (i in 2:length(p)) {
    aux[i] <- prod(
      zeta(alpha[i - 1], shift = p[i]) /
        zeta(alpha[i - 1], shift = p[i - 1]), aux[(i - 1)])
  }

  aux
}
