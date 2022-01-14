#' Read New Hampshire eXtended format for trees
#' @param fn Full path to the tree file.
#' @param txt If no file is specified, trees can also be passed as
#' a character scalar (see examples).
#' @return A list with the following elements:
#' - tree An object of class `ape`
#' - edge Edge annotations (length and other annotations)
#' - nhx A list of annotations NHX
#' @examples 
#' # Example directly extracted from
#' # https://sites.google.com/site/cmzmasek/home/software/forester/nhx
#' read_nhx(
#'   txt = "(((ADH2:0.1[&&NHX:S=human], ADH1:0.11[&&NHX:S=human]):0.05[&&NHX:S=primates:D=Y:B=100],
#'     ADHY:0.1[&&NHX:S=nematode],ADHX:0.12[&&NHX:S=insect]):0.1[&&NHX:S=metazoa:D=N],
#'     (ADH4:0.09[&&NHX:S=yeast],ADH3:0.13[&&NHX:S=yeast], ADH2:0.12[&&NHX:S=yeast],
#'     ADH1:0.11[&&NHX:S=yeast]):0.1 [&&NHX:S=Fungi])[&&NHX:D=N];"
#'     )
#' @export
#' @references "NHX - New Hampshire eXtended \[version 2.0\]", 
#' \url{https://en.wikipedia.org/wiki/Newick_format#New_Hampshire_X_format}
#' @family reading
read_nhx <- function(fn, txt) {
  
  if (!missing(fn) & !missing(txt))
    stop("Either of -fn- or -txt- should not be specified.", call. = FALSE)
  
  if (!missing(fn))
    txt <- paste(readLines(fn), collapse = "")
  
  
  # This pattern catches most of the information
  pattern <- "([(,)])\\s*([[:alnum:]_/-]*)\\s*([:][0-9]+[.]?[0-9]*)\\s*(\\[&&NHX[[:alnum:]_/:=\\s]+\\])?"
  
  # Capturing the patterns and splitting the data
  x <- gregexpr(txt, pattern = pattern, perl = TRUE)
  x <- regmatches(txt, x)
  
  # Creating a matrix with the data
  x <- regmatches(x[[1]], regexec(x[[1]], pattern = pattern, perl = TRUE))
  x <- do.call(rbind, x)
  
  # Do all have ids?
  noid <- which(x[,3] == "")
  if (length(noid)) {
    x[noid,3] <- sprintf("unnamed%04i", 1:length(noid))
  }
  
  for (i in noid) {
    txt <- sub(
      pattern     = x[i,1],
      replacement = paste0(x[i,2:4], collapse = ""),
      x           = txt,
      fixed       = TRUE
    )
  }
  
  # Is there any root?
  txt <- sub(
    pattern = "[)][;]$", replacement = ")root;", x = txt, perl = TRUE
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
    tree = ape::read.tree(text = txt),
    edge = dat[,-3L],
    nhx  = nhx
  )
  
}
