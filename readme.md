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
set.seed(123)
plot(sim_tree(10))
```

![](readme_files/figure-markdown_github/unnamed-chunk-1-1.png)

``` r
dat <- sim_annotated_tree(200, P=2)
as.apephylo(dat)
```

    ## 
    ## Phylogenetic tree with 200 tips and 199 internal nodes.
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
    ## [1] -253.0055
    ## 
    ## attr(,"class")
    ## [1] "phylo_LogLik"

MLE estimation
==============

``` r
# Using Artificial Bee Colony algorithm
ans0 <- phylo_mle(dat, method = "ABC")
ans0
```

    ## ESTIMATION OF ANNOTATED PHYLOGENETIC TREE
    ## ll: -196.0868,
    ## Method used: ABC (97 iterations)
    ## Leafs
    ##  # of Functions 2
    ##  # of 0:   100 (25%)
    ##  # of 1:   300 (75%)
    ## 
    ##          Estimate  Std. Error
    ##  psi[0]    0.0001      0.1070
    ##  psi[1]    0.1252      0.0323
    ##  mu[0]     0.1132      0.0697
    ##  mu[1]     0.0254      0.0129
    ##  Pi        0.9999      1.0527

``` r
# Plotting the path
plot(ans0)
```

![](readme_files/figure-markdown_github/MLE-1.png)
