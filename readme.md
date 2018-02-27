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
set.seed(198)
dat <- sim_annotated_tree(
  200, P=1, 
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
    ##      fun0000
    ## [1,]       1
    ## [2,]       0
    ## [3,]       0
    ## [4,]       1
    ## [5,]       1
    ## [6,]       1
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
    ## [6,]       0
    ## 
    ## ...(193 obs. omitted)...

Likelihood
----------

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
    ##  $ S : int [1:2, 1] 0 1
    ##  $ Pr: num [1:399, 1:2] 0.018 0.686 0.686 0.018 0.018 0.018 0.686 0.018 0.686 0.018 ...
    ##  $ ll: num -222
    ##  - attr(*, "class")= chr "phylo_LogLik"

Estimation
==========

``` r
# Using L-BFGS-B (MLE) to get an initial guess
ans0 <- aphylo_mle(dat)


# MCMC method
ans2 <- aphylo_mcmc(
  ans0$par, dat,
  prior = function(p) dbeta(p, 2,20),
  control = list(nbatch=1e4, burnin=100, thin=20, nchains=5))
ans2
```

    ## 
    ## ESTIMATION OF ANNOTATED PHYLOGENETIC TREE
    ## ll: -232.9435,
    ## Method used: mcmc (10000 iterations)
    ## Leafs
    ##  # of Functions 1
    ## 
    ##          Estimate  Std. Error
    ##  psi[0]    0.0956      0.0553
    ##  psi[1]    0.0549      0.0349
    ##  mu[0]     0.1212      0.0334
    ##  mu[1]     0.0917      0.0373
    ##  eta[0]    0.6414      0.0485
    ##  eta[1]    0.6954      0.0488
    ##  Pi        0.0983      0.0670

``` r
plot(ans2)
```

![](readme_files/figure-markdown_github/MLE-1.png)

``` r
# MCMC Diagnostics with coda
library(coda)
gelman.diag(ans2$hist)
```

    ## Potential scale reduction factors:
    ## 
    ##      Point est. Upper C.I.
    ## psi0       1.02       1.06
    ## psi1       1.01       1.01
    ## mu0        1.04       1.11
    ## mu1        1.01       1.03
    ## eta0       1.01       1.03
    ## eta1       1.00       1.01
    ## Pi         1.07       1.18
    ## 
    ## Multivariate psrf
    ## 
    ## 1.1

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
    ##         Mean      SD  Naive SE Time-series SE
    ## psi0 0.09556 0.05529 0.0011114       0.003289
    ## psi1 0.05488 0.03488 0.0007010       0.001396
    ## mu0  0.12119 0.03339 0.0006711       0.001381
    ## mu1  0.09165 0.03728 0.0007493       0.001598
    ## eta0 0.64137 0.04852 0.0009753       0.002561
    ## eta1 0.69541 0.04875 0.0009799       0.002465
    ## Pi   0.09827 0.06698 0.0013463       0.004882
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##          2.5%     25%     50%     75%  97.5%
    ## psi0 0.014520 0.05268 0.08818 0.12803 0.2251
    ## psi1 0.008309 0.02811 0.04762 0.07458 0.1390
    ## mu0  0.059180 0.09874 0.11978 0.14254 0.1905
    ## mu1  0.029218 0.06450 0.08839 0.11697 0.1721
    ## eta0 0.546001 0.60905 0.64115 0.67350 0.7345
    ## eta1 0.597325 0.66253 0.69601 0.72878 0.7861
    ## Pi   0.012759 0.04796 0.08324 0.13399 0.2684

Prediction
==========

``` r
pred <- prediction_score(ans2)
pred
```

    ## PREDICTION SCORE: ANNOTATED PHYLOGENETIC TREE
    ## Observed : 0.15 
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
