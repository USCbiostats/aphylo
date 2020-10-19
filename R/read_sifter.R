#' Read New Hampshire eXtended format for trees
#' @param fn Full path to the tree file
#' @return A list with the following elements:
#' - tree An object of class `ape`
#' - edge Edge annotations (length and other annotations)
#' - nhx A list of annotations NHX
read_nhx <- function(fn) {
  text <- paste(readLines(fn), collapse = "")
  
  # This pattern catches most of the information
  pattern <- "([(,)])([a-zA-Z0-9_]*)([:][0-9]+[.]?[0-9]*)(\\[&&NHX[a-zA-Z0-9_:=]+\\])?"
  
  # Capturing the patterns and splitting the data
  x <- gregexpr(text, pattern = pattern, perl = TRUE)
  x <- regmatches(text, x)
  
  # Creating a matrix with the data
  x <- regmatches(x[[1]], regexec(x[[1]], pattern = pattern, perl = TRUE))
  x <- do.call(rbind, x)
  
  # Do all have ids?
  noid <- which(x[,3] == "")
  if (length(noid)) {
    x[noid,3] <- sprintf("unnamed%04i", 1:length(noid))
  }
  
  for (i in noid) {
    text <- sub(
      pattern     = x[i,1],
      replacement = paste0(x[i,2:4], collapse = ""),
      x           = text,
      fixed       = TRUE
    )
  }
  
  # Is there any root?
  text <- sub(
    pattern = "[)][;]$", replacement = ")root;", x = text, perl = TRUE
  )
  
  dat <- x[,-c(1L, 2L)]
  
  # Capturing NHX fields
  nhx <- strsplit(dat[,3], split = "[:]|[=]")
  nhx <- tryCatch(lapply(nhx, function(n) {
    if (length(n) > 0) {
      n <- gsub(pattern = "(^\\[|\\]$)", replacement = "", x = n)
      n <- matrix(n[-1], ncol = 2, byrow = TRUE)
      structure(.Data = n[,2], names = n[,1])
    } else
      n
  }), error = function(e) e)
  
  if (inherits(nhx, "error")) 
    stop(
      "There was a problem when processing the &&NHS blocks.",
      " Possibly, not all the attributes have the right tag. Here is the error",
      ":\n", paste0(nhx, collapse=""), call. = FALSE
    )
  
  list(
    tree = ape::read.tree(text = text),
    edge = dat[,-3L],
    nhx  = nhx
  )
  
}

#' Read PLI files from SIFTER
#' @param fn Full path to the file
#' @param dropNAs Logical scalar. When `TRUE`, the function will discard any
#' protein that has no annotations.
#' @return A data table object including the following columns:
#' - name: Used to match UniProtKB data and GOA,
#' - number,
#' - go: A list of the GO annotations
#' - moc: Evidence code
#' - fam: Name of the family
#' @export
read_pli <- function(fn, dropNAs = TRUE) {

  ans <- xml2::as_list(xml2::read_html(fn))
  ans <- ans$html$body$family
  ans <- ans[which(names(ans) == "protein")]
  res <- lapply(ans, function(b) {
    
    # Extracting name
    name   <- unlist(b$proteinname)
    number <- unlist(b$proteinnumber)
    
    # And annotations
    go     <- unlist(b$gonumber)
    moc    <- unlist(b$moc)
    
    go <- if (is.null(go)) NA else unname(as.vector(go))
    moc <- if (is.null(moc)) NA else unname(as.vector(moc))
    
    # We are counting at least one annotation
    nann <- max(1, length(go))
    
    data.frame(
      name   = rep(unname(name), nann),
      number = rep(unname(number), nann),
      go     = go,
      moc    = moc,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
  })

  res        <- do.call(rbind, c(res, list(make.row.names = FALSE)))
  res$name   <- unlist(res$name)
  res$number <- unlist(res$number)
  res$go     <- strsplit(gsub("\\[|\\]|\\s", "", res$go), split = ",")
  res$moc    <- strsplit(gsub("\\[|\\]|\\s", "", res$moc), split = ",")
  # res[, c("name", "number") := list(unlist(name), unlist(number))]
  # res[, go := strsplit(gsub("\\[|\\]|\\s", "", go), split = ",")]
  # res[, moc := strsplit(gsub("\\[|\\]|\\s", "", moc), split = ",")]
  # 
  if (dropNAs)
    res <- res[!is.na(res$go),]
  
  go  <- res[["go"]]
  moc <- res[["moc"]]
  
  nrep <- sapply(go, length)
  nrep[nrep == 0] <- 1L
  
  res <- res[, which(!(colnames(res) %in% c("go", "moc")))]
  res <- res[rep(1:nrow(res), nrep),]
  
  cbind(
    res,
    data.frame(
      go  = unlist(go, recursive = TRUE),
      moc = unlist(moc, recursive = TRUE),
      stringsAsFactors = FALSE
      )
    )
  
  
}
