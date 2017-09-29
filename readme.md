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
    ##     6, 5, 4, 3, ...
    ## 
    ##   Internal nodes labels:
    ##     0, 2, 1, ...
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
    ## 
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
       offspring   = offspring, 
       psi = psi, mu = mu, Pi = pi_root)
)
```

    ## $ll
    ## [1] -57.35087
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
    ## ll:  -46.9586,
    ## Method used: L-BFGS-B (39 iterations)
    ## convergence: 0 (see ?optim)
    ## Leafs
    ##  # of Functions 1
    ##  # of 0:    25 (25%)
    ##  # of 1:    75 (75%)
    ## 
    ##          Estimate  Std. Error
    ##  psi[0]    0.2022      0.2499
    ##  psi[1]    0.0000      0.0609
    ##  mu[0]     0.0917      0.1366
    ##  mu[1]     0.0705      0.0382
    ##  Pi        1.0000      1.0702

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
    ## ll:  -42.5969,
    ## Method used: mcmc (10000 iterations)
    ## Leafs
    ##  # of Functions 1
    ##  # of 0:    25 (25%)
    ##  # of 1:    75 (75%)
    ## 
    ##          Estimate  Std. Error
    ##  psi[0]    0.1043      0.0613
    ##  psi[1]    0.0518      0.0309
    ##  mu[0]     0.1247      0.0562
    ##  mu[1]     0.0577      0.0219
    ##  Pi        0.1331      0.1029

``` r
# MCMC Diagnostics with coda
library(coda)
gelman.diag(ans2$hist)
```

    ## Potential scale reduction factors:
    ## 
    ##      Point est. Upper C.I.
    ## psi0       1.10       1.26
    ## psi1       1.01       1.02
    ## mu0        1.04       1.11
    ## mu1        1.01       1.03
    ## Pi         1.24       1.57
    ## 
    ## Multivariate psrf
    ## 
    ## 1.19

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
    ## psi0 0.10433 0.06128 0.0012319       0.004102
    ## psi1 0.05176 0.03086 0.0006203       0.001123
    ## mu0  0.12472 0.05623 0.0011303       0.003268
    ## mu1  0.05772 0.02189 0.0004400       0.000648
    ## Pi   0.13306 0.10285 0.0020674       0.009214
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##          2.5%     25%     50%     75%  97.5%
    ## psi0 0.013889 0.05625 0.09440 0.14451 0.2406
    ## psi1 0.007102 0.02800 0.04686 0.07014 0.1233
    ## mu0  0.030576 0.08445 0.12112 0.16152 0.2458
    ## mu1  0.021352 0.04225 0.05525 0.07099 0.1059
    ## Pi   0.017529 0.06576 0.10998 0.17293 0.4069

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
    ## Observed : 0.08 (73.72)
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
