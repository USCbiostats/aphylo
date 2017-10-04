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

``` r
# This datasets are included in the package
data("fakeexperiment")
data("faketree")

head(fakeexperiment)
```

    ##   LeafId f1 f2
    ## 1      3  0  0
    ## 2      4  0  1
    ## 3      5  1  0
    ## 4      6  1  1

``` r
head(faketree)
```

    ##      ParentId NodeId
    ## [1,]        1      3
    ## [2,]        1      4
    ## [3,]        2      5
    ## [4,]        2      6
    ## [5,]        0      1
    ## [6,]        0      2

``` r
O <- new_aphylo(
  annotations = fakeexperiment,
  edges       = faketree
)

O
```

    ## 
    ## A PARTIALLY ORDERED PHYLOGENETIC TREE
    ## 
    ##   # Internal nodes: 3
    ##   # Leaf nodes    : 4
    ## 
    ##   Leaf nodes labels: 
    ##     3, 4, 5, 6.
    ## 
    ##   Internal nodes labels:
    ##     0, 1, 2.
    ## 
    ## ANNOTATIONS:
    ##      f1 f2

``` r
as.apephylo(O)
```

    ## 
    ## Phylogenetic tree with 5 tips and 2 internal nodes.
    ## 
    ## Tip labels:
    ## [1] "3" "4" "5" "6" "0"
    ## Node labels:
    ## [1] "1" "2"
    ## 
    ## Rooted; includes branch lengths.

``` r
# We can visualize it
plot(O)
```

    ## Scale for 'fill' is already present. Adding another scale for 'fill',
    ## which will replace the existing scale.

![](readme_files/figure-markdown_github-ascii_identifiers/Get%20offspring-1.png)

``` r
plot_LogLike(O)
```

![](readme_files/figure-markdown_github-ascii_identifiers/Get%20offspring-2.png)

Simulating annoated trees
-------------------------

``` r
set.seed(1958)
dat <- sim_annotated_tree(
  100, P=1, 
  psi = c(0.05, 0.05),
  mu  = c(0.1, 0.05),
  Pi  = 1
  )

dat
```

    ## 
    ## A PARTIALLY ORDERED PHYLOGENETIC TREE
    ## 
    ##   # Internal nodes: 99
    ##   # Leaf nodes    : 100
    ## 
    ##   Leaf nodes labels: 
    ##     99, 100, 101, 102, 103, 104, ...
    ## 
    ##   Internal nodes labels:
    ##     0, 1, 2, 3, 4, 5, ...
    ## 
    ## ANNOTATIONS:
    ##      fun0000

Likelihood
----------

``` r
# Parameters and data
psi     <- c(0.020,0.010)
mu      <- c(0.04,.01)
pi_root <- .999

# Computing likelihood
with(dat, 
     LogLike(
       annotations = annotations, 
       offspring   = attr(edges, "offspring"), 
       psi = psi, mu = mu, Pi = pi_root)
)
```

    ## $ll
    ## [1] -68.26373
    ## 
    ## attr(,"class")
    ## [1] "phylo_LogLik"

Estimation
==========

``` r
# Using L-BFGS-B (MLE)
(ans0 <- aphylo_mle(dat))
```

    ## 
    ## ESTIMATION OF ANNOTATED PHYLOGENETIC TREE
    ## ll:  -50.1556,
    ## Method used: L-BFGS-B (30 iterations)
    ## convergence: 0 (see ?optim)
    ## Leafs
    ##  # of Functions 1
    ##  # of 0:    22 (22%)
    ##  # of 1:    78 (78%)
    ## 
    ##          Estimate  Std. Error
    ##  psi[0]    0.5394      0.2297
    ##  psi[1]    0.0724      0.0801
    ##  mu[0]     0.0369      0.1573
    ##  mu[1]     0.0677      0.0584
    ##  Pi        1.0000      1.0360

``` r
# Plotting loglike
plot_LogLike(ans0)
```

![](readme_files/figure-markdown_github-ascii_identifiers/MLE-1.png)

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
    ## ll:  -47.4017,
    ## Method used: mcmc (10000 iterations)
    ## Leafs
    ##  # of Functions 1
    ##  # of 0:    22 (22%)
    ##  # of 1:    78 (78%)
    ## 
    ##          Estimate  Std. Error
    ##  psi[0]    0.1340      0.0835
    ##  psi[1]    0.0831      0.0428
    ##  mu[0]     0.1623      0.0786
    ##  mu[1]     0.0578      0.0274
    ##  Pi        0.1312      0.1004

``` r
# MCMC Diagnostics with coda
library(coda)
gelman.diag(ans2$hist)
```

    ## Potential scale reduction factors:
    ## 
    ##      Point est. Upper C.I.
    ## psi0       1.11       1.27
    ## psi1       1.03       1.07
    ## mu0        1.06       1.16
    ## mu1        1.04       1.08
    ## Pi         1.11       1.26
    ## 
    ## Multivariate psrf
    ## 
    ## 1.11

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
    ## psi0 0.13401 0.08348 0.0016780       0.006577
    ## psi1 0.08315 0.04275 0.0008594       0.002059
    ## mu0  0.16226 0.07855 0.0015790       0.006112
    ## mu1  0.05780 0.02745 0.0005517       0.001168
    ## Pi   0.13120 0.10041 0.0020183       0.008871
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##         2.5%     25%     50%    75%  97.5%
    ## psi0 0.01861 0.07174 0.11816 0.1800 0.3394
    ## psi1 0.01533 0.05103 0.07904 0.1103 0.1775
    ## mu0  0.03132 0.10577 0.15593 0.2102 0.3364
    ## mu1  0.01539 0.03805 0.05467 0.0737 0.1191
    ## Pi   0.01810 0.06455 0.10694 0.1690 0.4257

``` r
plot(ans2$hist)
```

![](readme_files/figure-markdown_github-ascii_identifiers/MCMC-1.png)![](readme_files/figure-markdown_github-ascii_identifiers/MCMC-2.png)

Prediction
==========

``` r
pred <- prediction_score(ans2)
pred
```

    ## PREDICTION SCORE: ANNOTATED PHYLOGENETIC TREE
    ## Observed : 0.08 (81.86)
    ## Random   : 0.25 (243.26)
    ## Best     : 0.00 (0.00)
    ## Worse    : 1.00 (973.05)
    ## ---------------------------------------------------------------------------
    ## Values between 0 and 1, 0 been best. Absolute scores in parenthesis.

``` r
plot(pred)
```

![](readme_files/figure-markdown_github-ascii_identifiers/Predict-1.png)

Misc
====

During the development process, we decided to allow the user to choose what 'tree-reader' function he would use, in particular, between using either the rncl R package or ape. For such we created a short benchmark that compares both functions [here](playground/ape_now_supports_singletons.md).
