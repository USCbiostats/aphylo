#' Functional balance of a tree
#' 
#' This function computes the distance between .5 and the observed proportion
#' of ones for each function in a tree.
#' 
#' @param phy An object of class [aphylo] or [multiAphylo]
#' 
#' @details 
#' Functional balance is defined as follows
#' 
#' \deqn{%
#' P^{-1}\sum_{p}\left(1 - \left|0.5 - N^{-1}\sum_n a_{np}\right|\right)
#' }{%
#' mean(1 - abs(.5 - colMeans(A)))
#' }
#'  
#' Where `A` is the matrix of annotations.
#'  
#' With values ranging between 0 and 1, one been perfect balance, this is, equal
#' number of zeros and ones in the annotations. In the case of multiple functions,
#' as noted in the formula, the balance is the average across functions. 
#' 
#' @return If `phy` is an object of class `phylo`, a single scalar, otherwise,
#' it returns a vector of length [`Ntrees`]`(phy)`.
#' 
#' @examples 
#' x <- raphylo(20, P = 2)
#' balance_ann(x)
#' 
#' balance_ann(c(x, x))
#' 
#' @export
balance_ann <- function(phy) UseMethod("balance_ann")

#' @export
balance_ann.aphylo <- function(phy) {
  
  A <- phy$tip.annotation
  A[A==9L] <- NA_integer_
  
  mean((.5 - abs(.5 - colMeans(A, na.rm = TRUE)))/.5)
  
}

#' @export
balance_ann.multiAphylo <- function(phy) {
  
  sapply(phy, balance_ann)
  
}

#' @export
balance_ann.default <- function(phy) {
  
  stop(
    "No defined method to calculate functional balance in objects of class '",
    paste(class(phy), collapse = ", "), "'", call. = FALSE
    )
  
}
