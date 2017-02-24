# scoping
rm(list = ls())

# D and f are defined in the Global env
D <- "I'm global"

f <- function() print(D)

# The D in this function is defined in the environment of G
g <- function(fun0) {
  D <- "I'm local"
  
  # This should obviate the local D
  fun0()
  
  # Not this, as fun1 is defined in the same env than 
  fun1 <- function() print(D) 
  fun1()
}

# This should print:
# [1] "I'm global"
# [1] "I'm local"
g(f)
