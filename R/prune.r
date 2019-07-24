#' @export
#' @rdname new_aphylo_pruner
#' @param ... Further arguments passed to the method
new_aphylo_pruner <- function(...) UseMethod("new_aphylo_pruner")

#' @export
#' @param annotation a list of length `N` (annotations).
#' @rdname new_aphylo_pruner
new_aphylo_pruner.matrix <- function(edgelist, annotation, Ntype = NULL, ...) {
  
  n <- max(edgelist)
  P <- ncol(annotation)
  
  if (ncol(edgelist) != 2L)
    stop("`edgelist` must be a two-column matrix.", call. = FALSE)
  
  if (nrow(annotation) != n)
    stop("`nrow(annotation) != nrow(edge`", call. = FALSE)
  
  if (is.null(Ntype))
    Ntype <- rep(1L, n)
  else if (length(Ntype) != n)
    stop("`nrow(annotation) != nrow(edge`", call. = FALSE)
    
  # Checking the complete cases
  nannotated <- annotation
  nannotated[nannotated == 9L] <- NA
  nannotated <- sum(complete.cases(nannotated))
  
  # Splitting the data
  annotation <- lapply(seq_len(n), function(i) as.vector(annotation[i, ]))
  
  new_aphylo_pruner.(
    edgelist = list(edgelist[, 1L] - 1L, edgelist[, 2L] - 1L),
    A        = annotation,
    Ntype    = Ntype
    )
  
}

#' @export
#' @rdname new_aphylo_pruner
#' @param x An object of class [aphylo].
new_aphylo_pruner.aphylo <- function(x, ...) {
  
  annotation <- with(x, rbind(tip.annotation, node.annotation))
  annotation <- lapply(seq_len(nrow(annotation)), function(i) annotation[i, ])
  
  new_aphylo_pruner.(
    edgelist = list(x$tree$edge[, 1L] - 1L, x$tree$edge[, 2L] - 1L),
    A        = annotation,
    Ntype    = rep(1L, length(annotation)), 
    nannotated = x$Ntips.annotated
  )
  
}

#' @export
#' @rdname new_aphylo_pruner
new_aphylo_pruner.multiAphylo <- function(x, ...) {
  
  structure(
    lapply(x, new_aphylo_pruner, ...),
    class = "multiAphylo_pruner"
  )
  
}