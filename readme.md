Phylogenetic tree
================

Install
-------

Using devtools

``` r
devtools::install_github("gvegayon/phylogenetic")
```

Example
-------

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
O <- get_offsprings(
    experiment, "LeafId", 
    tree, "NodeId", "ParentId"
)

# There is no nice print method for now
lapply(O, head)
```

    ## $offsprings
    ## $offsprings[[1]]
    ## integer(0)
    ## 
    ## $offsprings[[2]]
    ## integer(0)
    ## 
    ## $offsprings[[3]]
    ## integer(0)
    ## 
    ## $offsprings[[4]]
    ## integer(0)
    ## 
    ## $offsprings[[5]]
    ## integer(0)
    ## 
    ## $offsprings[[6]]
    ## integer(0)
    ## 
    ## 
    ## $noffsprings
    ## [1] 0 0 0 0 0 0
    ## 
    ## $edgelist
    ##      NodeId ParentId
    ## [1,]      2        0
    ## [2,]      6        2
    ## [3,]      7        6
    ## [4,]      9        7
    ## [5,]     13        9
    ## [6,]     16       13

``` r
# We can visualize it
plot(O, vertex.size=5, vertex.label=NA)
```

![](readme_files/figure-markdown_github/Get%20offsprings-1.png)

Likelihood
----------

``` r
# Parameters and data
Z <- as.matrix(experiment[,-4])

psi     <- c(0.020,0.010)
mu      <- c(0.004,.001)
pi_root <- c(1-0.1,.1)

# Computing likelihood
ll1   <- LogLike(Z, O$offsprings, O$noffsprings, psi, mu, pi_root)
lapply(ll1, head)
```

    ## $S
    ##      [,1] [,2] [,3]
    ## [1,]    0    0    0
    ## [2,]    1    0    0
    ## [3,]    0    1    0
    ## [4,]    1    1    0
    ## [5,]    0    0    1
    ## [6,]    1    0    1
    ## 
    ## $PI
    ##       [,1]
    ## [1,] 0.729
    ## [2,] 0.081
    ## [3,] 0.081
    ## [4,] 0.009
    ## [5,] 0.081
    ## [6,] 0.009
    ## 
    ## $PSI
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
    ## [1,] 0.01 0.01 0.01 0.01 0.98 0.98 0.98 0.98
    ## [2,] 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
    ## [3,] 0.01 0.01 0.98 0.98 0.01 0.01 0.98 0.98
    ## [4,] 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
    ## [5,] 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
    ## [6,] 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
    ## 
    ## $Pr
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
    ## [1,] 0.01 0.01 0.01 0.01 0.98 0.98 0.98 0.98
    ## [2,] 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
    ## [3,] 0.01 0.01 0.98 0.98 0.01 0.01 0.98 0.98
    ## [4,] 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
    ## [5,] 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
    ## [6,] 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
    ## 
    ## $ll
    ## [1] -47.39684
