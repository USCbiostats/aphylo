#' Indexing aphylo objects
#' @param x An object of class [aphylo].
#' @param i,j Integer vector. Indices of genes or functions.
#' @param value Integer vector. Replacing values, can be either `c(0, 1, 9, NA)`.
#' @param drop Logical scalar. When `TRUE`, the function returns a matrix of 
#' annotations. Otherwise an object of class `aphylo`.
#' @details The subsetting method allows selecting one or more annotations from
#' the [aphylo] object. Whenever `i` is specified, then aphylo returns the corresponding
#' annotations.
#' @name aphylo-index
#' @returns 
#' - When indexing with `i`: A data frame with the annotations of the
#'   selected genes.
#' - When only indexing with `j` (`drop = FALSE`): An `aphylo` object with the selected sets of
#'   annotations.
#' - When only indexing with `j` (`drop = TRUE`): A data.frame with the selected
#'   annotations.
#' - When indexing on both `i` and `j`: A data.frame with the selected genes and annotations.
#' @examples 
#' set.seed(12312)
#' atree <- raphylo(50, P = 4)
#' atree[1:10,]
#' atree[,2:3]
#' atree[, 2:3, drop = TRUE]
#' atree[1:10, 2:3]
#' 
#' 
NULL

#' @export
#' @rdname aphylo-index
`[.aphylo` <- function(x, i, j, drop = FALSE) {
  
  if (missing(i) & missing(j)) {
    
    stop("You should either specify i or j when indexing aphylo objects.", call. = FALSE)
    
  } else if (missing(i)) {
    
    if (drop) {
      return(rbind(x$tip.annotation, x$node.annotation)[, j, drop=FALSE])
    } 
    
    x$tip.annotation  <- x$tip.annotation[, j, drop=FALSE]
    x$node.annotation <- x$node.annotation[, j, drop=FALSE]
    
    # Updating the pruning sequence
    x$reduced_pseq <- reduce_pseq(
      pseq      = x$pseq,
      A         = with(x, rbind(tip.annotation, node.annotation)),
      offspring = list_offspring(x)
      )
    
    return(x)
    
  } else if (missing(j)) {
    
    if (drop)
      warning("drop = TRUE has no effect when subsetting genes.", call. = FALSE)
    
    return(rbind(x$tip.annotation, x$node.annotation)[i, , drop = FALSE])
    
  } else {
    
    return(rbind(x$tip.annotation, x$node.annotation)[i, j, drop = FALSE])
    
  }
  
  
}

#' @export
#' @rdname aphylo-index
`[<-.aphylo` <- function(x, i, j, value) {
  
  # Checking all values
  value[is.na(value)] <- 9L
  test <- which(!(value %in% c(0L, 1L, 9L)))
  if (length(test)) 
    stop(
      "Replacing annotations must be either c(0, 1, 9, NA).",
      " The following values are not valid: ",
      paste(value[test], collapse=", "), ".", call. = FALSE
      )
  
  if (missing(j))
    j <- 1L:Nann(x)
  if (missing(i))
    i <- 1L:Nnode(x, internal.only = FALSE)
  
  # Sorting
  nt <- Ntip(x)
  tips_ids  <- which(i <= nt)
  nodes_ids <- which(i > nt)
  
  if (length(tips_ids))
    x$tip.annotation[i[tips_ids], j] <- value[tips_ids]
  
  if (length(nodes_ids))
    x$node.annotation[i[nodes_ids] - nt, j] <- value[nodes_ids]
  
  x
}

