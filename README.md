# pwpldis Package

<!-- badges: start -->
[![R-CMD-check](https://github.com/njerezlillo/pwpldis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/njerezlillo/pwpldis/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
![Lifecycle: experimental](https://img.shields.io/badge/Lifecycle-Experimental-orange)
<!-- badges: end -->

This package provides tools for fitting the piecewise discrete power-law model, a flexible statistical framework for modeling data that exhibits power-law behavior in different segments. It is applicable to fields such as complex networks, empirical distributions, and heavy-tailed phenomena. The package includes:

- Base Functions: Implementation of various functions related to the model, such as density, survival, hazard, random number generation, and more.
- Maximum Likelihood Estimation: A method for estimating model parameters, including change points and scaling parameters.
- Bootstrap Bias Correction: A resampling procedure for correcting estimation bias and quantifying uncertainty.

## Progress status

- [x] Set up the package structure  
- [x] Write functions  
- [x] Document functions
- [x] Write examples for each function in the package
- [x] Check the documentation
- [x] Publish on GitHub  
- [ ] Complete the "Example" section on GitHub
- [ ] Distribute on CRAN

## Installation

You can install the package using:

``` r
# install.packages("devtools")
devtools::install_github("njerezlillo/pwpldis")
```

## Example

Load the required package `pwpldis` which contains functions for fitting and 
simulating discrete piecewise power-law distributions.

``` r
library(pwpldis)
```

Setting a fixed random seed ensures that the results of the simulation are reproducible. 
This is particularly useful for simulations, model validation, and debugging.

``` r
set.seed(2025)
```

Define the breakpoints `p` that separate different power-law regimes in the dataset. 
The `alpha` values correspond to the exponents governing each segment of the power-law model.

``` r
p <- c(1, 10)
alpha <- c(1.5, 3.5)
```

The function `rpwpldis` generates a random sample from a discrete piecewise power-law distribution. Here, we generate a dataset of 1000 observations based on the specified 
breakpoints and exponents.

``` r
df <- rpwpldis(1000, p, alpha)
```

The function `fit_pwpldis` estimates the parameters of a discrete piecewise power-law model. 
The argument `nbreak = 1` indicates that the model should be fitted with exactly one breakpoint.
The function returns an estimated change point location and fitted power-law exponents.

``` r
fit_1 <- fit_pwpldis(df, nbreak = 1)
```

The bootstrap procedure is used to estimate the variability of the model parameters. 
Resampling is performed with replacement, and parameters are re-estimated for each resampled dataset. The predefined breakpoint is taken from the fitted model `fit_1$tau_1` to ensure consistency.

``` r
boot_1 <- boot_pwpldis(df, brks = fit_1$tau_1)
```

Bias correction is performed using a standard bootstrap adjustment method. 
The bias is estimated by comparing the original parameter estimates to the bootstrap 
averages. The corrected estimates are obtained as: 
`corrected_estimate = 2 * original_estimate - mean(bootstrap_estimates)`.

``` r
bias_corrected_params <- 2 * fit_1[, 3:4] - apply(boot_1[, 3:4], 2, mean)
```

The ECDF of the simulated dataset is plotted to visually assess the distribution of the data.
This function plots the empirical distribution function using black points for clarity.

``` r
plot(ecdf(df), cex = 0.5, main = "Empirical vs. Fitted CDF", xlab = "x", ylab = "CDF")
```

To compare the empirical distribution with the fitted theoretical model, we compute 
the theoretical cumulative distribution function (CDF). The function `ppwpldis` evaluates 
the CDF based on the estimated parameters.

``` r
### Overlay the Theoretical Cumulative Distribution Function (CDF)
cdf <- Vectorize(function(x) ppwpldis(x, c(1, 10), c(1.48, 3.49)))

# Overlay the theoretical CDF on the empirical plot using red points.
points(unique(df), cdf(unique(df)), col = "red", cex = 0.5)

# Add a legend to distinguish between the empirical and theoretical distributions.
legend("bottomright", legend = c("Empirical CDF", "Fitted CDF"), col = c("black", "red"), 
pch = c(1, 1))
```

## Citation

To cite `pwpldis` package in publications, please use the following format:

Jerez-Lillo N (2025). *pwpldis: Piecewise Discrete Power-Law Model*. R package version 1.0.0, [https://github.com/njerezlillo/pwpldis](https://github.com/njerezlillo/pwpldis).

For LaTeX users, the corresponding BibTeX entry is:

```bibtex
@Manual{
  title = {pwpldis: Piecewise Discrete Power-Law Model},
  author = {Nixon Jerez-Lillo},
  year = {2025},
  note = {R package version 1.0.0},
  url = {https://github.com/njerezlillo/pwpldis},
}
```
