Phylogenetic tree
================

[![Travis-CI Build Status](https://travis-ci.org/USCbiostats/phylogenetic.svg?branch=master)](https://travis-ci.org/USCbiostats/phylogenetic) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/USCbiostats/phylogenetic?branch=master&svg=true)](https://ci.appveyor.com/project/USCbiostats/phylogenetic) [![Coverage Status](https://img.shields.io/codecov/c/github/USCbiostats/phylogenetic/master.svg)](https://codecov.io/github/USCbiostats/phylogenetic?branch=master)

Install
-------

Using devtools

``` r
devtools::install_github("USCbiostats/phylogenetic")
```

Reading data
------------

``` r
library(phylogenetic)
```

``` r
# This datasets are included in the package
data("experiment")
data("tree")

head(experiment)
```

    ##   f01 f02 f03 LeafId
    ## 1   9   9   1     16
    ## 2   9   9   9     20
    ## 3   9   1   9     21
    ## 4   9   9   9     34
    ## 5   9   9   9     39
    ## 6   9   9   9     50

``` r
head(tree)
```

    ##   NodeId TypeId ParentId
    ## 1      2      0        0
    ## 2      6      1        2
    ## 3      7      0        6
    ## 4      9      0        7
    ## 5     13      0        9
    ## 6     16      2       13

``` r
O <- get_offspring(
    experiment, "LeafId", 
    tree, "NodeId", "ParentId"
)

# There is no nice print method for now
as.phylo(O)
```

    ## 
    ## Phylogenetic tree with 175 tips and 152 internal nodes.
    ## 
    ## Tip labels:
    ##  leaf001, leaf002, leaf003, leaf004, leaf005, leaf006, ...
    ## 
    ## Unrooted; includes branch lengths.

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
```

Likelihood
----------

``` r
# Parameters and data

psi     <- c(0.020,0.010)
mu      <- c(0.04,.01)
pi_root <- c(1-0.5,.5)

# Computing likelihood
LogLike(dat$experiment, dat$offspring, dat$noffspring, psi, mu, pi_root)
```

    ## $ll
    ## [1] -282.5007
    ## 
    ## attr(,"class")
    ## [1] "phylo_LogLik"

MLE estimation
==============

``` r
# Using Artificial Bee Colony algorithm
ans0 <- phylo_mle(rep(.2,5), dat, useABC = TRUE, control = list(maxCycle = 50))

# Plotting the path
plot(ans0)
```

![](readme_files/figure-markdown_github/MLE-1.png)
