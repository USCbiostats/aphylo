library(testthat)
library(aphylo)

Sys.setenv("R_TESTS" = "")
test_check("aphylo")
