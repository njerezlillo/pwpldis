% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cj_pwpldis.R
\name{Cj_pwpldis}
\alias{Cj_pwpldis}
\title{Normalizing Constants for the Discrete Piecewise Power-Law Model}
\usage{
Cj_pwpldis(p, alpha)
}
\arguments{
\item{p}{A numeric vector specifying the minimum value \eqn{\tau_{(0)}} and
change points \eqn{\tau_{(1)}, \ldots, \tau_{(k)}} that define the partitions,
i.e., \eqn{p = (\tau_{(0)}, \tau_{(1)}, \ldots, \tau_{(k)})}. It must be sorted
in increasing order.}

\item{alpha}{A numeric vector of scaling parameters
\eqn{\alpha_{(1)}, \ldots, \alpha_{(k + 1)}} corresponding to each partition.
All values must be greater than 1.}
}
\value{
A numeric vector containing the normalizing constants for each
partition of the model. By definition, the first value returned
is 1, and the last is 0.
}
\description{
This function computes the normalizing constants required for the discrete
piecewise power-law model, ensuring that the probability mass function sums
to one over the defined partitions.
}
\details{
The constants are given by:
\deqn{C_0=1, \quad C_{i}=\displaystyle\prod_{h=1}^{i} \frac{\zeta(\alpha_h,\tau_{(h)})}{\zeta(\alpha_{h},\tau_{(h-1)})},}
with \eqn{i=1,\ldots,k.}

The function assumes that \code{p} and \code{alpha} have the same length.

The calculation involves the Hurwitz zeta function from the VGAM package.
}
\examples{
# Example with increasing parameters
Cj_pwpldis(c(1, 3), c(1.5, 2.0))

# Example with decreasing parameters
Cj_pwpldis(c(1, 3, 5), c(2.5, 2.0, 1.5))

# Example with non-equidistant partitions
Cj_pwpldis(c(1, 2, 4, 7), c(2.0, 2.5, 1.5, 2.0))

}
\references{
Jerez-Lillo, N., Rodrigues, F. A., Ferreira, P. H., & Ramos, P. L. (2025).
Beyond the Power Law: Estimation, Goodness-of-Fit, and a Semiparametric
Extension in Complex Networks. arXiv preprint arXiv:2311.11200. Available at:
\url{https://arxiv.org/abs/2311.11200}
}
