#' Create an `aphylo` object with partial annotations
#' 
#' @param tree An object of class `phylo`.
#' @param annotations A [data.frame] with annotations. The first column should
#' be the gene id (see details).
#' @param types A [data.frame] with types. Just like the annotations, the first
#' column should be the gene id.
#' 
#' @details Each row in the the `annotations` data frame passed to this function
#' must have a unique row per gene, and one column per function (GO term). The id
#' of each gene must match the labels in the `tree` object. Missing genes are
#' annotated with `NA` (9).
#' 
#' In the case of `types`, while tips can also be annotated with a type, which
#' should be either 0, duplication, or 1, speciation, only internal nodes
#' are required. Tip types are ignored.
#' 
#' @return An object of class [aphylo].
#' 
#' @examples 
#' # Generating a test dataset
#' set.seed(1371)
#' x <- raphylo(20)
#' 
#' # Extracting the tree and annotations
#' tree <- x$tree
#' 
#' anno <- with(x, rbind(tip.annotation, node.annotation))
#' anno <- data.frame(id = with(tree, c(tip.label, node.label)), anno)
#' 
#' types <- data.frame(id = tree$node.label, x$node.type)
#' 
#' # Creating a aphylo tree without node types
#' aphylo_from_data_frame(tree, anno)
#' 
#' # Now including types
#' aphylo_from_data_frame(tree, anno, types)
#' 
#' # Dropping some data
#' aphylo_from_data_frame(tree, anno[sample.int(nrow(anno), 10),])
#' 
#' @export
#' @family Data management functions
aphylo_from_data_frame <- function(tree, annotations, types = NULL) {
  
  # Checking the types
  if (!inherits(tree, "phylo"))
    stop(
      "-tree- must be an object of class \"phylo\".",
      " The object passed is of class \"", class(tree), "\".", call. = FALSE
      )
  
  if (!inherits(annotations, "data.frame"))
    stop(
      "-annotations- must be an object of class \"data.frame\".",
      " The object passed is of class \"", class(annotations), "\".", call. = FALSE
    )
  
  if (!is.null(types) && !inherits(types, "data.frame"))
    stop(
      "-types- must be an object of class \"data.frame\".",
      " The object passed is of class \"", class(types), "\".", call. = FALSE
    )
  
  if (ncol(annotations) < 2L)
    stop(
      "The -annotations- data frame must have at least two columns, one for the",
      " gene label (id) and another for the corresponding annotations. The ",
      "-annotations- data frame has only ", ncol(annotations), " columns.",
      call. = FALSE
      )
  
  if (!is.null(types) && ncol(types) != 2L)
    stop(
      "The -types- data frame must have two columns, one for the",
      " gene label (id) and another for the corresponding type. The ",
      "-types- data frame has ", ncol(types), " columns.",
      call. = FALSE
    )
  
  # Checking annotations, first, are there any labels?
  tip_labels  <- tree$tip.label
  node_labels <- tree$node.label
  
  # Basic checks of the data
  if (length(c(tip_labels, node_labels)) == 0)
    stop("Neither the tips nor the nodes have labels in this tree.",
         call. = FALSE)
  
  if (!all(annotations[,1] %in% c(tip_labels, node_labels)))
    stop(
      "Not all the entries in -annotations- are present in the labels of",
      " the phylo object (tree).", call. = FALSE
    )
  
  # Preparing tip annotations --------------------------------------------------
  if (is.null(tip_labels))
    stop("This tree has no tip.labels.", call. = FALSE)
    
  # Which entries are tip annotations
  tip_idx <- match(tip_labels, as.character(annotations[, 1L]))
  
  # Preparing the annotation matrix
  tip_annotations <- matrix(
    NA_integer_, ncol = ncol(annotations) - 1, nrow = length(tip_labels),
    dimnames = list(NULL, colnames(annotations)[-1L])
    )
  
  not_miss <- which(!is.na(tip_idx))
  
  tip_annotations[not_miss,] <- as.matrix(annotations[tip_idx[not_miss], -1L])
  
  # Working on node annotations ------------------------------------------------
  if (!is.null(node_labels)) {
    
    # Which entries are tip annotations
    node_idx <- match(node_labels, as.character(annotations[, 1L]))
    
    # Preparing the annotation matrix
    node_annotation <- matrix(
      NA_integer_, ncol = ncol(annotations) - 1, nrow = length(node_labels),
      dimnames = list(NULL, colnames(annotations)[-1L])
    )
    
    not_miss <- which(!is.na(node_idx))
    
    node_annotation[not_miss,] <- as.matrix(annotations[node_idx[not_miss], -1L])
    
  } else
    node_annotation <- NULL
  
  # Working on types ------------------------------------------------------------
  if (!is.null(types)) {
    
    if (is.null(node_labels))
      stop(
        "When specifying -types-, the tree's nodes must be labeled. This tree",
        " has no labeled nodes.",
        call. = FALSE
        )
    
    if (!all(types[,1] %in% c(tip_labels, node_labels)))
      stop(
        "Not all the entries in -types- are present in the labels of",
        " the phylo object (tree).", call. = FALSE
      )
    
    # Node types (all should be specified)
    test <- which(!(node_labels %in% types[,1]))
    if (length(test))
      stop(
        "When specifying -types-, all node labels should be present in the",
        " first column. The following are missing:\n",
        node_labels[test]
        )
    
    # Which entries are tip annotations
    node_idx <- match(node_labels, as.character(types[, 1L]))
    
    # Preparing the annotation matrix
    node_types <- matrix(
      NA_integer_, ncol = 1L, nrow = length(node_labels),
      dimnames = list(NULL, "type")
    )
    
    not_miss <- which(!is.na(node_idx))
    
    node_types[not_miss,] <- types[node_idx[not_miss], -1L]
   
  } else
    node_types <- NULL
  
  new_aphylo(
    tree            = tree,
    tip.annotation  = tip_annotations,
    node.annotation = node_annotation,
    tip.type        = NULL,
    node.type       = node_types
  )
  
}

