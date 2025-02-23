% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot_pwpldis.R
\name{boot_pwpldis}
\alias{boot_pwpldis}
\title{Bootstrap for Piecewise Power-Law Distribution}
\usage{
boot_pwpldis(
  time,
  nsim = 100,
  brks = NULL,
  nbreak = 1,
  exclude_int = NULL,
  min_pt_tail = 5,
  max_set = 1000,
  tol = 1e-04,
  parallel = FALSE,
  mc.core = 4,
  ...
)
}
\arguments{
\item{time}{A numeric vector of observations.}

\item{nsim}{An integer specifying the number of bootstrap simulations.
Default is 100.}

\item{brks}{A numeric vector of change points for the piecewise power-law model.
If \code{NULL}, the change points are estimated.}

\item{nbreak}{An integer specifying the number of change points in the model.
Default is 1.}

\item{exclude_int}{A numeric vector of two values specifying an interval where
change points should not be placed (e.g., \code{exclude_int = c(5, Inf)} excludes
change points greater than 5).}

\item{min_pt_tail}{An integer specifying the minimum number of observations
required in the last segment to ensure stable parameter estimation. Default is 5.}

\item{max_set}{An integer specifying the maximum number of possible change point
combinations considered for model fitting. Default is 1000.}

\item{tol}{A numeric value specifying the tolerance for combining close change points.
Default is 1e-4.}

\item{parallel}{A logical value indicating whether to use parallel processing
for bootstrap simulations. Default is \code{FALSE}.}

\item{mc.core}{An integer specifying the number of cores to use in parallel processing.
Default is 4.}

\item{...}{Additional arguments passed to the underlying \code{fit_pwpldis} function.}
}
\value{
A \code{data.frame} with the results of the bootstrap procedure, including the
following columns:
\describe{
\item{taus}{Estimated change points from each bootstrap iteration.}
\item{alphas}{Estimated scaling parameters for each segment from each
bootstrap iteration.}
\item{likelihood}{The log-likelihood in each bootstrap iteration.}
\item{AIC}{The Akaike Information Criterion in each bootstrap iteration.}
\item{BIC}{The Bayesian Information Criterion in each bootstrap iteration.}
}
}
\description{
This function performs a bootstrap procedure to estimate the distribution of
estimators for a piecewise power-law model. It resamples the dataset multiple
times and fits the model parameters for each bootstrap sample.
}
\details{
This function performs a bootstrap procedure by resampling the original dataset
(\code{time}) with replacement for \code{nsim} iterations. For each resample, the function
fits a piecewise power-law model using either the specified or estimated change points,
and the model parameters are recorded.

The \code{exclude_int} parameter prevents change points from being placed in a specific interval.
The \code{min_pt_tail} parameter ensures that the last segment has a sufficient number of
observations for stable parameter estimation.

After the bootstrap simulations, the results are aggregated and returned,
including the estimated change points and scaling parameters for each segment.

This function was inspired by the implementation of the \code{PWEXP}
package in \code{R}. Visit their GitHub repository
\href{https://github.com/zjph602xtc/PWEXP}{PWEXP}.
}
\examples{
# Initial specification
set.seed(2025)
p <- c(1, 3, 5)
alpha <- c(1.5, 2, 3)

# Generate a dataset from the discrete piecewise power-law model
df <- rpwpldis(300, p, alpha)

# Fit the piecewise power-law distribution with a predefined change point at 3
fit_1 <- fit_pwpldis(df, breakpoint = 3)

# Perform the bootstrap procedure with a predefined change point at 3
boot_1 <- boot_pwpldis(df, brks = fit_1$tau_1)

# Calculate the bias-corrected parameters
2 * fit_1[, 3:4] - apply(boot_1[, 3:4], 2, mean)

# Calculate the 95\% confidence intervals for the bootstrap estimates of the model parameters
apply(boot_1[, 3:4], 2, function(x) quantile(x, c(0.025, 0.975)))

# Fit the model with 1 change point estimated from the data
fit_2 <- fit_pwpldis(df, nbreak = 1)

# Perform the bootstrap procedure with 1 change point estimated from the data
boot_2 <- boot_pwpldis(df, nbreak = 1)

# Calculate the bias-corrected parameters for the second model
2 * fit_2[, 3:4] - apply(boot_2[, 3:4], 2, mean)

# Calculate the 95\% confidence intervals for the bootstrap estimates of the second model
apply(boot_2[, 3:4], 2, function(x) quantile(x, c(0.025, 0.975)))

}
\references{
Jerez-Lillo, N., Rodrigues, F. A., Ferreira, P. H., & Ramos, P. L. (2025).
Beyond the Power Law: Estimation, Goodness-of-Fit, and a Semiparametric
Extension in Complex Networks. arXiv preprint arXiv:2311.11200. Available at:
\url{https://arxiv.org/abs/2311.11200}
}
\seealso{
\link{fit_pwpldis}
}
