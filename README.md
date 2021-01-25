
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CopSens

<!-- badges: start -->
<!-- badges: end -->

`CopSens` implements a copula-based sensitivity analysis method for
unobserved confounding in multi-treatment inference.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JiajingZ/CopSens")
```

The dependency `pcaMethods` from [Bioconductor](http://bioconductor.org)
may fail to automatically install, which would result in a warning
similar to:

    ERROR: dependency ‘pcaMethods’ is not available for package ‘CopSens’

Then, please first install the `pcaMethods` manually by

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pcaMethods")
```

## Basic Usage

``` r
# load package #
library(CopSens)
#> Loading required package: tidyverse
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──
#> ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
#> ✓ tibble  3.0.5     ✓ dplyr   1.0.3
#> ✓ tidyr   1.1.2     ✓ stringr 1.4.0
#> ✓ readr   1.4.0     ✓ forcats 0.5.0
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
```

### Example for analysis with Gaussian outcomes

``` r
# load data #
y <- GaussianT_GaussianY$y
tr <- subset(GaussianT_GaussianY, select = -c(y))

# execute worst-case calibration #
est_df1 <- gcalibrate(y = y, tr = tr, t1 = tr[1,], t2 = tr[2,],
                      calitype = "worstcase", R2 = c(0.6, 1))
#> Fitting the latent confounder model by PPCA with default.
#> 1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:
#> Observed outcome model fitted by simple linear regression with default.
#> Worst-case calibration executed.
# visualize #
plot_estimates(est_df1)
```

<img src="man/figures/README-gaussian-outcome-example-1.png" width="85%" />

``` r
# execute multivariate calibration #
est_df2 <- gcalibrate(y = y, tr = tr, t1 = tr[1:10,], t2 = tr[11:20,],
                      calitype = "multicali")
#> Fitting the latent confounder model by PPCA with default.
#> 1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:
#> Observed outcome model fitted by simple linear regression with default.
#> Multivariate calibration executed.
# visualize #
plot_estimates(est_df2)
```

<img src="man/figures/README-gaussian-outcome-example-2.png" width="85%" />

``` r
# execute user-specified calibration #
est_df3 <- gcalibrate(y = y, tr = tr, t1 = tr[1:2,], t2 = tr[3:4,],
                      calitype = "null", gamma = c(0.96, -0.29, 0),
                      R2 = c(0.3, 0.7, 1))
#> Fitting the latent confounder model by PPCA with default.
#> 1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:
#> Observed outcome model fitted by simple linear regression with default.
#> User-specified calibration executed.
# visualize #
plot_estimates(est_df3)
```

<img src="man/figures/README-gaussian-outcome-example-3.png" width="85%" />

### Example for analysis with binary outcomes

``` r
# load data #
y <- GaussianT_BinaryY$y
tr <- subset(GaussianT_BinaryY, select = -c(y))
t1 <- tr[1:5,]
t2 <- rep(0, times = ncol(tr))

# calibrate #
est_df <- bcalibrate(y = y, tr = tr, t = rbind(t1, t2),
                     gamma = c(1.27, -0.28, 0), R2 = c(0.5, 0.7))
#> Fitting the latent confounder model by PPCA with default.
#> 1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:1:2:3:4:5:6:7:8:9:10:
#> Observed outcome model fitted by simple probit model with default.
#> R2 =  0.5 , calibrating observation 1  2  3  4  5  6  
#> R2 =  0.7 , calibrating observation 1  2  3  4  5  6
# calculate risk ratio estimator #
rr_df <- est_df[1:5,] / as.numeric(est_df[6,])
# visualize #
plot_estimates(rr_df)
```

<img src="man/figures/README-binary-outcome-example-1.png" width="85%" />
