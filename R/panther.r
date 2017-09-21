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
#' read_panther(path)
#' 
#' @name panther-tree
#' @aliases PANTHER PANTHERDB
NULL

#' @export
#' @param tree.reader Function that will be used to read the tree file.
#' It can be either \code{ape::read.tree} or \code{rncl::read_newick_phylo}.
#' @rdname panther-tree
read_panther <- function(x, tree.reader = ape::read.tree, ...) {
  # Reading the data-in
  x  <- readLines(x)
  
  # Obtaining extra info and processing internal nodes labels
  rgxp <- "(?:[:])?([0-9.]+)?\\[\\&\\&NHX:Ev=([0-9><]{3})(?::S=([a-zA-Z_.-]+))?:ID=([a-zA-Z0-9]+)\\]"
  
  dat <- stringr::str_match_all(x, rgxp)[[1]]
  
  # Rewriting the file so that labels of inner nodes can be read in
  for (i in 1:nrow(dat))
    x[1] <- stringr::str_replace(
      string      = x[1],
      pattern     = stringr::fixed(dat[i,1]),
      replacement = ifelse(is.na(dat[i,2]), dat[i,5], paste(dat[i,5],dat[i,2], sep=":"))
      )
  
  # Getting the labels
  labs <- data.frame(
    id       = stringr::str_extract(x[-1], "^.+(?=\\:)"),
    label    = stringr::str_extract(x[-1], "(?<=\\:).+(?=\\;$)"),
    stringsAsFactors = FALSE
  )
  
  # Reading the tree
  # tree <- ape::read.tree(text=x[1], ...)
  tmptree <- tempfile()
  write(x[1], tmptree)
  tree <- tree.reader(tmptree, ...)
  file.remove(tmptree)
  
  tree$tip.label <- paste(
    tree$tip.label,
    labs$label[match(tree$tip.label, labs$id)],
    sep=":"
  )
  
  # Creating a nice data-frame
  ans <- data.frame(
    branch_length    = as.numeric(dat[,2]),
    type             = ifelse(dat[,3] == "0>1", "S",
                              ifelse(dat[,3] == "1>0", "D", "T")),
    ancestor         = dat[,4],
    row.names        = dat[,5], 
    stringsAsFactors = FALSE
  )
  
  # Which ones are duplication nodes
  ans$duplication <- ifelse(ans$type %in% c("D", "T"), TRUE, FALSE)
  
  # Sorting and returning
  list(
    tree = tree,
    internal_nodes_annotations  = ans[order(as.integer(gsub("[a-zA-Z]+","",rownames(ans)))),]
  )
}

#' @rdname panther-tree
#' @export
read.panther <- read_panther