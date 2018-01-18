aphylo: Statistical Inference of Annotated Phylogenetic Trees
================

[![Travis-CI Build Status](https://travis-ci.org/USCbiostats/aphylo.svg?branch=master)](https://travis-ci.org/USCbiostats/aphylo) [![Build status](https://ci.appveyor.com/api/projects/status/1xpgpv10yifojyab?svg=true)](https://ci.appveyor.com/project/gvegayon/phylogenetic) [![Coverage Status](https://codecov.io/gh/USCbiostats/aphylo/branch/master/graph/badge.svg)](https://codecov.io/gh/USCbiostats/aphylo)

The `aphylo` R package implements estimation and data imputation methods for Functional Annotations in Phylogenetic Trees. The core function consists on the computation of the log-likelihood of observing a given phylogenetic tree with functional annotation on its leafs, and probabilities associated to gain and loss of functionalities, including probabilities of experimental misclassification. Furthermore, the log-likelihood is computed using peeling algorithms, which required developing and implementing efficient algorithms for re-coding and preparing phylogenetic tree data so that can be used with the package. Finally, `aphylo` works smoothly with popular tools for analysis of phylogenetic data such as `ape` R package, "Analyses of Phylogenetics and Evolution".

The package is under MIT License, and is been developed by the Computing and Software Cores of the Biostatistics Division's NIH Project Grant (P01) at the Department of Preventive Medicine at the University of Southern California.

Install
-------

This package depends on another on-development R package, the [`amcmc`](https://github.com/USCbiostats/amcmc). So first you need to install it:

``` r
devtools::install_github("USCbiostats/amcmc")
```

Then you can install the `aphylo` package

``` r
devtools::install_github("USCbiostats/aphylo")
```

Reading data
------------

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

    ## Scale for 'fill' is already present. Adding another scale for 'fill',
    ## which will replace the existing scale.

![](readme_files/figure-markdown_github/Get%20offspring-1.png)

``` r
plot_LogLike(O)
```

![](readme_files/figure-markdown_github/Get%20offspring-2.png)

Simulating annoated trees
-------------------------

``` r
set.seed(1958)
dat <- sim_annotated_tree(
  200, P=1, 
  psi = c(0.05, 0.05),
  mu  = c(0.1, 0.1),
  Pi  = .5
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
    ##      fun0000
    ## [1,]       0
    ## [2,]       0
    ## [3,]       1
    ## [4,]       0
    ## [5,]       1
    ## [6,]       0
    ## 
    ## ...(194 obs. omitted)...
    ## 
    ## 
    ##  Internal node annotations:
    ##      fun0000
    ## [1,]       1
    ## [2,]       1
    ## [3,]       1
    ## [4,]       1
    ## [5,]       1
    ## [6,]       1
    ## 
    ## ...(193 obs. omitted)...

Likelihood
----------

``` r
# Parameters and data
psi     <- c(0.020,0.010)
mu      <- c(0.04,.01)
pi_root <- .999

# Computing likelihood
str(LogLike(dat, psi = psi, mu = mu, Pi = pi_root))
```

    ## List of 3
    ##  $ S : int [1:2, 1] 0 1
    ##  $ Pr: num [1:399, 1:2] 0.98 0.98 0.02 0.98 0.02 0.98 0.02 0.02 0.02 0.02 ...
    ##  $ ll: num -158
    ##  - attr(*, "class")= chr "phylo_LogLik"

Estimation
==========

``` r
# Using L-BFGS-B (MLE)
(ans0 <- aphylo_mle(dat))
```

    ## Warning in sqrt(diag(x$varcovar)): NaNs produced

    ## 
    ## ESTIMATION OF ANNOTATED PHYLOGENETIC TREE
    ## ll: -135.5278,
    ## Method used: L-BFGS-B (41 iterations)
    ## convergence: 0 (see ?optim)
    ## Leafs
    ##  # of Functions 1
    ## 
    ##          Estimate  Std. Error
    ##  psi[0]    1.0000      0.2089
    ##  psi[1]    1.0000      0.5945
    ##  mu[0]     0.7675      0.5698
    ##  mu[1]     0.6075      0.4436
    ##  Pi        0.0000         NaN

``` r
# Plotting loglike
plot_LogLike(ans0)
```

![](readme_files/figure-markdown_github/MLE-1.png)

``` r
# MCMC method
ans2 <- aphylo_mcmc(
  ans0$par, dat,
  prior = function(p) dbeta(p, 2,20),
  control = list(nbatch=1e4, burnin=100, thin=20, nchains=5))
ans2
```

    ## 
    ## ESTIMATION OF ANNOTATED PHYLOGENETIC TREE
    ## ll: -115.5646,
    ## Method used: mcmc (10000 iterations)
    ## Leafs
    ##  # of Functions 1
    ## 
    ##          Estimate  Std. Error
    ##  psi[0]    0.1098      0.0863
    ##  psi[1]    0.0831      0.1121
    ##  mu[0]     0.0893      0.0909
    ##  mu[1]     0.1496      0.0669
    ##  Pi        0.1123      0.0708

``` r
# MCMC Diagnostics with coda
library(coda)
gelman.diag(ans2$hist)
```

    ## Potential scale reduction factors:
    ## 
    ##      Point est. Upper C.I.
    ## psi0       1.03       1.08
    ## psi1       1.03       1.06
    ## mu0        1.02       1.05
    ## mu1        1.01       1.02
    ## Pi         1.14       1.34
    ## 
    ## Multivariate psrf
    ## 
    ## 1.12

``` r
summary(ans2$hist)
```

    ## 
    ## Iterations = 120:10000
    ## Thinning interval = 20 
    ## Number of chains = 5 
    ## Sample size per chain = 495 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##         Mean      SD Naive SE Time-series SE
    ## psi0 0.10983 0.08633 0.001735       0.007050
    ## psi1 0.08308 0.11212 0.002254       0.012421
    ## mu0  0.08926 0.09088 0.001827       0.009196
    ## mu1  0.14961 0.06687 0.001344       0.005693
    ## Pi   0.11227 0.07078 0.001423       0.005244
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##          2.5%     25%     50%     75%  97.5%
    ## psi0 0.020432 0.06493 0.09309 0.12996 0.3584
    ## psi1 0.008046 0.03252 0.05610 0.08608 0.5097
    ## mu0  0.017525 0.04982 0.07065 0.09546 0.4537
    ## mu1  0.079970 0.11770 0.13905 0.16248 0.3496
    ## Pi   0.013754 0.05811 0.09981 0.15649 0.2744

``` r
plot(ans2$hist)
```

![](readme_files/figure-markdown_github/MCMC-1.png)![](readme_files/figure-markdown_github/MCMC-2.png)

Prediction
==========

``` r
pred <- prediction_score(ans2)
pred
```

    ## PREDICTION SCORE: ANNOTATED PHYLOGENETIC TREE
    ## Observed : 0.14 
    ## Random   : 0.25 
    ## ---------------------------------------------------------------------------
    ## Values standarized to range between 0 and 1, 0 being best.

``` r
plot(pred)
```

![](readme_files/figure-markdown_github/Predict-1.png)

Misc
====

During the development process, we decided to allow the user to choose what 'tree-reader' function he would use, in particular, between using either the rncl R package or ape. For such we created a short benchmark that compares both functions [here](playground/ape_now_supports_singletons.md).
