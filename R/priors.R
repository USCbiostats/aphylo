#' Default priors for [aphylo_mcmc]
#' 
#' Convenient wrappers to be used with the `aphylo` estimation methods.
#' 
#' @param shape1,shape2,... Arguments passed to [stats::dbeta]
#' @return In the case of `bprior`, a wrapper of the function [stats::dbeta].
#' `uprior` returns a function `function(p) 1` (the uniform prior)
#' @examples 
#' bprior(1, 9)
#' uprior()
#' @export
bprior <- function(shape1=1, shape2=9, ...) {
  
  function(p) {
    
    stats::dbeta(p, shape1 = shape1, shape2 = shape2, ...)
    
  }
  
}

#' @rdname bprior
#' @export
uprior <- function() {
  function(p) 1
}

