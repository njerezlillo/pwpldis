
# pwpldis Package

<!-- badges: start -->
[![R-CMD-check](https://github.com/njerezlillo/pwpldis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/njerezlillo/pwpldis/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
![Lifecycle: experimental](https://img.shields.io/badge/Lifecycle-Experimental-orange)
<!-- badges: end -->

This package provides a set of tools for fitting the piecewise discrete power-law model, a flexible statistical framework for modeling data that follows a power-law behavior in different segments. It includes:

- Maximum Likelihood Estimation: Methods to estimate model parameters, including breakpoints and power-law exponents.
- Bootstrap Bias Correction: A resampling procedure to correct estimation bias and quantify uncertainty.

The package is particularly useful for analyzing datasets where the power-law behavior changes at specific points, making it applicable to fields such as complex networks, empirical distributions, and heavy-tailed phenomena.

## Progress status

- [x] Set up the package structure  
- [x] Write functions  
- [x] Document functions
- [ ] Write examples for each function in the package
- [ ] Check the documentation
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

This is a basic example which shows you how to solve a common problem:

``` r
library(pwpldis)
## basic example code
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
