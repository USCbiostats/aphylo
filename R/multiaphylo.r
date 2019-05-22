#' Building Lists of Annotated Trees
#' 
#' This is equivalent to what [ape::c.phylo] does.
#' 
#' @param ... One or several object of class `aphylo` or `multiAPhylo`.
#' @name multiAphylo
#' @examples 
#' data(fakeexperiment)
#' data(faketree)
#' ans <- new_aphylo(fakeexperiment[,2:3], faketree)
#' c(ans, ans)
NULL

#' @export
#' @rdname multiAphylo
c.aphylo <- function(...) {
  
  ans  <- list(...)
  
  if (length(ans) == 1)
    return(ans)
  
  # Checking classes. We don't allow combining objects other than
  # aphylo or multiAphylo.
  are_aphylo      <- sapply(ans, is.aphylo)
  are_multiAphylo <- sapply(ans, is.multiAphylo)
  
  if (any(!are_aphylo & !are_multiAphylo))
    stop(
      "Some elements of `...` are neither `ahylo` nor `multiAphylo` ",
      "objects.", call. = FALSE
    )
  
  # Concatenating
  are_multiAphylo <- which(are_multiAphylo)
  are_aphylo      <- which(are_aphylo)
  
  structure(
    c(
      if (length(are_multiAphylo)) {
        unclass(do.call(c, ans[[are_multiAphylo]]))
      } else NULL,
      if (length(are_aphylo)) {
        unclass(ans[are_aphylo])
      } else NULL
    ),
    class = "multiAphylo"
  )
  
}

#' @export
#' @rdname multiAphylo
c.multiAphylo <- c.aphylo

print.multiAphylo <- function(x, ...) {
  
  N <- length(x)
  cat(N, "annotated phylogenetic", ifelse(N > 1, "trees\n", "tree\n"))
  
  invisible(x)
  
}
