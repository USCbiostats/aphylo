aphylo: Statistical Inference of Annotated Phylogenetic Trees
================

[![Travis-CI Build Status](https://travis-ci.org/USCbiostats/aphylo.svg?branch=master)](https://travis-ci.org/USCbiostats/aphylo) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/USCbiostats/aphylo?branch=master&svg=true)](https://ci.appveyor.com/project/USCbiostats/aphylo) [![Coverage Status](https://img.shields.io/codecov/c/github/USCbiostats/aphylo/master.svg)](https://codecov.io/github/USCbiostats/aphylo?branch=master)

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

    ##   f1 f2 LeafId
    ## 1  0  0      3
    ## 2  0  1      4
    ## 3  1  0      5
    ## 4  1  1      6

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
  leafidvar   = "LeafId",
  edges       = faketree
)

# There is no nice print method for now
as.apephylo(O)
```

    ## 
    ## Phylogenetic tree with 4 tips and 3 internal nodes.
    ## 
    ## Tip labels:
    ## [1] "leaf001" "leaf002" "leaf003" "leaf004"
    ## 
    ## Rooted; includes branch lengths.

``` r
# We can visualize it
plot(O)
```

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
  100, P=1, 
  psi = c(0.05, 0.05),
  mu  = c(0.1, 0.05),
  Pi  = 1
  )

as.apephylo(dat)
```

    ## 
    ## Phylogenetic tree with 100 tips and 99 internal nodes.
    ## 
    ## Tip labels:
    ##  leaf001, leaf002, leaf003, leaf004, leaf005, leaf006, ...
    ## 
    ## Rooted; includes branch lengths.

Likelihood
----------

``` r
# Parameters and data
psi     <- c(0.020,0.010)
mu      <- c(0.04,.01)
pi_root <- c(1-0.5,.5)

# Computing likelihood
with(dat, 
     LogLike(
       annotations = annotations, 
       offspring   = offspring, 
       noffspring  = noffspring, 
       psi = psi, mu = mu, Pi = pi_root)
)
```

    ## $ll
    ## [1] -58.03654
    ## 
    ## attr(,"class")
    ## [1] "phylo_LogLik"

Estimation
==========

``` r
# Using L-BFGS-B (MLE)
(ans0 <- phylo_mle(dat))
```

    ## ESTIMATION OF ANNOTATED PHYLOGENETIC TREE
    ## ll:  -46.9587,
    ## Method used: L-BFGS-B (31 iterations)
    ## convergence: 0 (see ?optim)
    ## Leafs
    ##  # of Functions 1
    ##  # of 0:    25 (25%)
    ##  # of 1:    75 (75%)
    ## 
    ##          Estimate  Std. Error
    ##  psi[0]    0.2022      0.2499
    ##  psi[1]    0.0001      0.0609
    ##  mu[0]     0.0917      0.1366
    ##  mu[1]     0.0705      0.0382
    ##  Pi        0.9999      1.0699

``` r
# Using ABC (MLE)
(ans1 <- phylo_mle(dat, method="ABC"))
```

    ## ESTIMATION OF ANNOTATED PHYLOGENETIC TREE
    ## ll:   47.0033,
    ## Method used: ABC (55 iterations)
    ## Leafs
    ##  # of Functions 1
    ##  # of 0:    25 (25%)
    ##  # of 1:    75 (75%)
    ## 
    ##          Estimate  Std. Error
    ##  psi[0]    0.2582      0.2434
    ##  psi[1]    0.0001      0.0573
    ##  mu[0]     0.0527      0.1335
    ##  mu[1]     0.0711      0.0379
    ##  Pi        0.9999      1.0263

``` r
# MCMC method
ans2 <- phylo_mcmc(ans0$par, dat, control = list(nbatch=1e4, burnin=100, thin=20,
                                         nchains=10 # Will run 10 chains
                                         ))
ans2
```

    ## ESTIMATION OF ANNOTATED PHYLOGENETIC TREE
    ## ll:  -49.3773,
    ## Method used: mcmc (10000 iterations)
    ## Leafs
    ##  # of Functions 1
    ##  # of 0:    25 (25%)
    ##  # of 1:    75 (75%)
    ## 
    ##          Estimate  Std. Error
    ##  psi[0]    0.2247      0.1424
    ##  psi[1]    0.0512      0.0423
    ##  mu[0]     0.1365      0.0880
    ##  mu[1]     0.0703      0.0367
    ##  Pi        0.6135      0.2743

``` r
# MCMC Diagnostics with coda
library(coda)
gelman.diag(ans2$hist)
```

    ## Potential scale reduction factors:
    ## 
    ##      Point est. Upper C.I.
    ## psi0       1.13       1.25
    ## psi1       1.00       1.01
    ## mu0        1.07       1.14
    ## mu1        1.02       1.04
    ## Pi         1.25       1.49
    ## 
    ## Multivariate psrf
    ## 
    ## 1.3

``` r
summary(ans2$hist)
```

    ## 
    ## Iterations = 120:20000
    ## Thinning interval = 20 
    ## Number of chains = 10 
    ## Sample size per chain = 995 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##         Mean      SD  Naive SE Time-series SE
    ## psi0 0.22470 0.14241 0.0014277       0.009918
    ## psi1 0.05123 0.04234 0.0004244       0.001054
    ## mu0  0.13655 0.08797 0.0008819       0.004386
    ## mu1  0.07033 0.03669 0.0003679       0.001733
    ## Pi   0.61347 0.27432 0.0027501       0.034867
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##          2.5%     25%     50%     75%  97.5%
    ## psi0 0.012907 0.11077 0.21055 0.31956 0.5376
    ## psi1 0.001841 0.01924 0.04146 0.07175 0.1606
    ## mu0  0.007645 0.06562 0.12584 0.19439 0.3300
    ## mu1  0.020567 0.04678 0.06524 0.08674 0.1473
    ## Pi   0.052025 0.40843 0.65396 0.85133 0.9871

``` r
plot(ans2$hist)
```

![](readme_files/figure-markdown_github/MLE-1.png)![](readme_files/figure-markdown_github/MLE-2.png)
