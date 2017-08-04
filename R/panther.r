#' Reads PANTHER db trees 
#' 
#' The PANTHER Project handles a modified version of newick tree files which,
#' besides of the tree structure, includes the type of node and ancestor
#' labels. This function is a wrapper of \code{\link[ape:read.tree]{read.tree}}.
#' 
#' @param x Character scalar. Full path to the panther file.
#' @param ... Further arguments passed to \code{\link[ape:read.tree]{read.tree}}.
#' @return
#' 
#' A list consisting of a data.frame and a \code{phylo} object. The
#' data.frame has the following columns:
#'  
#' \item{branch_length}{Numeric vector. Length of the branch to its parent node.}
#' \item{type}{Character vector. Can be either \code{"S"} (speciation), \code{"D"}
#' (duplication), or \code{"T"} (horizontal transfer).}
#' \item{ancestor}{Character vector. Name of the ancestor.}
#' 
#' The nodeids can be identified using the rownames.
#' 
#' @examples
#' path <- system.file("tree.tree", package="aphylo")
#' read_panther_tree(path)
#' 
#' @name panther-tree
#' @export
read_panther_tree <- function(x, ...) {
  # Reading the data-in
  x  <- readLines(x, n = 1L)
  
  tree <- ape::read.tree(text=x, ...)
  
  # Matching the expression
  rgxp <- "(?:[:])?([0-9.]+)?\\[\\&\\&NHX:Ev=([0-9><]{3})(?::S=([a-zA-Z_.-]+))?:ID=([a-zA-Z0-9]+)\\]"
  dat <- stringr::str_match_all(x, rgxp)[[1]][,-1]
  
  # Creating a nice data-frame
  ans <- data.frame(
    branch_length    = as.numeric(dat[,1]),
    type             = ifelse(dat[,2] == "0>1", "S",
                              ifelse(dat[,2] == "1>0", "D", "T")),
    ancestor         = dat[,3],
    row.names        = dat[,4],
    stringsAsFactors = FALSE
  )
  
  # Which ones are duplication nodes
  ans$duplication <- ifelse(ans$type %in% c("D", "T"), TRUE, FALSE)
  
  # Sorting and returning
  list(
    tree = tree,
    dat  = ans[order(as.integer(gsub("[a-zA-Z]+","",rownames(ans)))),]
  )
}