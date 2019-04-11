aphylo: Statistical Inference of Annotated Phylogenetic Trees
================

[![](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Travis-CI Build
Status](https://travis-ci.org/USCbiostats/aphylo.svg?branch=master)](https://travis-ci.org/USCbiostats/aphylo)
[![Build
status](https://ci.appveyor.com/api/projects/status/1xpgpv10yifojyab?svg=true)](https://ci.appveyor.com/project/gvegayon/phylogenetic)
[![Coverage
Status](https://codecov.io/gh/USCbiostats/aphylo/branch/master/graph/badge.svg)](https://codecov.io/gh/USCbiostats/aphylo)

The `aphylo` R package implements estimation and data imputation methods
for Functional Annotations in Phylogenetic Trees. The core function
consists on the computation of the log-likelihood of observing a given
phylogenetic tree with functional annotation on its leafs, and
probabilities associated to gain and loss of functionalities, including
probabilities of experimental misclassification. Furthermore, the
log-likelihood is computed using peeling algorithms, which required
developing and implementing efficient algorithms for re-coding and
preparing phylogenetic tree data so that can be used with the package.
Finally, `aphylo` works smoothly with popular tools for analysis of
phylogenetic data such as `ape` R package, “Analyses of Phylogenetics
and Evolution”.

The package is under MIT License, and is been developed by the Computing
and Software Cores of the Biostatistics Division’s NIH Project Grant
(P01) at the Department of Preventive Medicine at the University of
Southern California.

## Install

This package depends on another on-development R package, the
[`amcmc`](https://github.com/USCbiostats/amcmc). So first you need to
install it:

``` r
devtools::install_github("USCbiostats/amcmc")
```

Then you can install the `aphylo` package

``` r
devtools::install_github("USCbiostats/aphylo")
```

## Reading data

``` r
library(aphylo)
```

    ## Loading required package: ape

``` r
# This datasets are included in the package
data("fakeexperiment")
data("faketree")

head(fakeexperiment)
```

    ##      LeafId f1 f2
    ## [1,]      1  0  0
    ## [2,]      2  0  1
    ## [3,]      3  1  0
    ## [4,]      4  1  1

``` r
head(faketree)
```

    ##      ParentId NodeId
    ## [1,]        6      1
    ## [2,]        6      2
    ## [3,]        7      3
    ## [4,]        7      4
    ## [5,]        5      6
    ## [6,]        5      7

``` r
O <- new_aphylo(
  tip.annotation = fakeexperiment[,2:3],
  tree           = faketree
)

O
```

    ## 
    ## Phylogenetic tree with 4 tips and 3 internal nodes.
    ## 
    ## Tip labels:
    ## [1] 1 2 3 4
    ## Node labels:
    ## [1] 5 6 7
    ## 
    ## Rooted; no branch lengths.
    ## 
    ##  Tip (leafs) annotations:
    ##   f1 f2
    ## 1  0  0
    ## 2  0  1
    ## 3  1  0
    ## 4  1  1
    ## 
    ##  Internal node annotations:
    ##   f1 f2
    ## 5  9  9
    ## 6  9  9
    ## 7  9  9

``` r
as.phylo(O)
```

    ## 
    ## Phylogenetic tree with 4 tips and 3 internal nodes.
    ## 
    ## Tip labels:
    ## [1] 1 2 3 4
    ## Node labels:
    ## [1] 5 6 7
    ## 
    ## Rooted; no branch lengths.

``` r
# We can visualize it
plot(O)
```

![](README_files/figure-gfm/Get%20offspring-1.png)<!-- -->

``` r
plot_logLik(O)
```

    ## No parameters were specified. Default will be used instead.

![](README_files/figure-gfm/Get%20offspring-2.png)<!-- -->

## Simulating annoated trees

``` r
set.seed(198)
dat <- raphylo(
  200, P=2, 
  psi = c(0.05, 0.05),
  mu  = c(0.1, 0.1),
  eta = c(.7, .95),
  Pi  = .4
  )

dat
```

    ## 
    ## Phylogenetic tree with 200 tips and 199 internal nodes.
    ## 
    ## Tip labels:
    ##  1, 2, 3, 4, 5, 6, ...
    ## Node labels:
    ##  201, 202, 203, 204, 205, 206, ...
    ## 
    ## Rooted; no branch lengths.
    ## 
    ##  Tip (leafs) annotations:
    ##      fun0000 fun0001
    ## [1,]       1       1
    ## [2,]       0       1
    ## [3,]       0       9
    ## [4,]       1       1
    ## [5,]       1       1
    ## [6,]       1       1
    ## 
    ## ...(194 obs. omitted)...
    ## 
    ## 
    ##  Internal node annotations:
    ##      fun0000 fun0001
    ## [1,]       1       0
    ## [2,]       1       0
    ## [3,]       1       1
    ## [4,]       1       1
    ## [5,]       1       1
    ## [6,]       0       1
    ## 
    ## ...(193 obs. omitted)...

## Likelihood

``` r
# Parameters and data
psi     <- c(0.020,0.010)
mu      <- c(0.04,.01)
eta     <- c(.7, .9)
pi_root <- .05

# Computing likelihood
str(LogLike(dat, psi = psi, mu = mu, eta = eta, Pi = pi_root))
```

    ## List of 3
    ##  $ S : int [1:4, 1:2] 0 1 0 1 0 0 1 1
    ##  $ Pr: num [1:399, 1:4] 0.000324 0.012348 0.203056 0.000324 0.000324 ...
    ##  $ ll: num -399
    ##  - attr(*, "class")= chr "phylo_LogLik"

# Estimation

``` r
# Using L-BFGS-B (MLE) to get an initial guess
ans0 <- aphylo_mle(dat ~ psi + mu + Pi + eta)
```

    ## No parameters were specified. Default will be used instead.

``` r
# MCMC method
ans2 <- aphylo_mcmc(
  dat ~ mu + psi + eta,
  prior = function(p) dbeta(p, 2,20),
  control = list(nsteps=1e4, burnin=100, thin=20, nchains=5))
```

    ## No parameters were specified. Default will be used instead.

    ## Warning: A single initial point has been passed via `initial`: c(0.1, 0.1,
    ## 0.1, 0.1, 0.9, 0.9). The values will be recycled.

    ## Convergence has been reached with 1100 steps (50 final count of observations).

``` r
ans2
```

    ## 
    ## ESTIMATION OF ANNOTATED PHYLOGENETIC TREE
    ## 
    ##  Call: aphylo_mcmc(model = dat ~ mu + psi + eta, priors = function(p) dbeta(p, 
    ##     2, 20), control = list(nsteps = 10000, burnin = 100, thin = 20, 
    ##     nchains = 5))
    ##  LogLik (unnormalized): -425.5128
    ##  Method used: mcmc (1100 steps)
    ##  # of Leafs: 200
    ##  # of Functions 2
    ##          Estimate  Std. Err.
    ##  psi0    0.0879    0.0475
    ##  psi1    0.0438    0.0272
    ##  mu0     0.1139    0.0262
    ##  mu1     0.0929    0.0218
    ##  eta0    0.6770    0.0349
    ##  eta1    0.8006    0.0320

``` r
plot(ans2)
```

![](README_files/figure-gfm/MLE-1.png)<!-- -->

``` r
# MCMC Diagnostics with coda
library(coda)
gelman.diag(ans2$hist)
```

    ## Potential scale reduction factors:
    ## 
    ##      Point est. Upper C.I.
    ## psi0       1.03       1.09
    ## psi1       1.02       1.08
    ## mu0        1.03       1.10
    ## mu1        1.01       1.05
    ## eta0       1.08       1.22
    ## eta1       1.01       1.04
    ## 
    ## Multivariate psrf
    ## 
    ## 1.11

``` r
summary(ans2$hist)
```

    ## 
    ## Iterations = 120:1100
    ## Thinning interval = 20 
    ## Number of chains = 5 
    ## Sample size per chain = 50 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##         Mean      SD Naive SE Time-series SE
    ## psi0 0.08789 0.04745 0.003001       0.006267
    ## psi1 0.04376 0.02715 0.001717       0.002065
    ## mu0  0.11386 0.02619 0.001656       0.003464
    ## mu1  0.09293 0.02183 0.001381       0.001936
    ## eta0 0.67703 0.03495 0.002210       0.003062
    ## eta1 0.80062 0.03205 0.002027       0.002872
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##          2.5%     25%     50%     75%  97.5%
    ## psi0 0.015681 0.05758 0.07854 0.11911 0.1945
    ## psi1 0.003866 0.02378 0.03856 0.06099 0.1095
    ## mu0  0.064390 0.09719 0.11336 0.13285 0.1652
    ## mu1  0.054532 0.07710 0.09142 0.10847 0.1320
    ## eta0 0.611915 0.65315 0.67380 0.70114 0.7434
    ## eta1 0.737275 0.78161 0.80163 0.82253 0.8594

# Prediction

``` r
pred <- prediction_score(ans2)
pred
```

    ## PREDICTION SCORE: ANNOTATED PHYLOGENETIC TREE
    ## Observed : 0.01 
    ## Random   : 0.42 
    ## ---------------------------------------------------------------------------
    ## Values scaled to range between 0 and 1, 0 being best.

``` r
plot(pred)
```

![](README_files/figure-gfm/Predict-1.png)<!-- -->![](README_files/figure-gfm/Predict-2.png)<!-- -->

# Misc

During the development process, we decided to allow the user to choose
what ‘tree-reader’ function he would use, in particular, between using
either the rncl R package or ape. For such we created a short benchmark
that compares both functions
[here](playground/ape_now_supports_singletons.md).
