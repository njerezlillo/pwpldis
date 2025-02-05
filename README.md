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
