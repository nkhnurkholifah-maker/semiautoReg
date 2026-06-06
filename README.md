
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semiautoReg

<!-- badges: start -->

<!-- badges: end -->

The goal of `semiautoReg` is to provide a simple workflow for
semiparametric regression with autocorrelated errors.

This package combines spline-based flexible regression and
autocorrelation correction in one R package structure. It is developed
as part of the Statistical Programming course.

## Main Features

- Automatic identification of linear and nonlinear predictors
- Automatic AR(q) order selection
- Group SCAD penalized estimation
- B-spline based nonlinear modeling
- Diagnostic plotting tools
- Prediction for new observations

  
## Installation

The package can be installed using:

``` r
# install.packages("remotes")
remotes::install_github("nkhnurkholifah-maker/semiautoReg")
```

## Example

```r
library(semiautoReg)

data(macro_turkey)

prep <- semiauto_prepare(
  data = macro_turkey,
  y = "inflasi",
  x = c("kurs","m3","industri","kredit","brent"),
  time = "Tarih"
)

fit <- semiauto_fit(prepared = prep)

summary(fit)
```
