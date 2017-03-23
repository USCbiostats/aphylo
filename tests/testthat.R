library(testthat)
library(phylogenetic)

Sys.setenv("R_TESTS" = "")
test_check("phylogenetic")
