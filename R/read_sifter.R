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
#' @family reading
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

#' Write pli files used by SIFTER
#' @param protein_name,protein_number,go_number,moc Vectors of the same length
#' @param family_id Character scalar. Name of the family
#' @param file Character scalar passed to [cat].
#' @export
#' @returns 
#' A string with the XML file.
#' @examples 
#' set.seed(882)
#' atree <- raphylo(5)
#' write_pli(
#'   family_id      = "a family",
#'   protein_name   = atree$tree$tip.label,
#'   protein_number = 1:Ntip(atree),
#'   go_number      = "GO:123123123123"
#' )
#' # Possible outcome:
#' #<?xml version="1.0"?>
#' #<Family>
#' #  <FamilyID>a family</FamilyID>
#' #  <Protein>
#' #    <ProteinName>1</ProteinName>
#' #    <ProteinNumber>1</ProteinNumber>
#' #    <GONumber>[GO:123123123123]</GONumber>
#' #    <MOC>[EXP]</MOC>
#' #  </Protein>
#' #  <Protein>
#' #    <ProteinName>2</ProteinName>
#' #    <ProteinNumber>2</ProteinNumber>
#' #    <GONumber>[GO:123123123123]</GONumber>
#' #    <MOC>[EXP]</MOC>
#' #  </Protein>
#' #  <Protein>
#' #    <ProteinName>3</ProteinName>
#' #    <ProteinNumber>3</ProteinNumber>
#' #    <GONumber>[GO:123123123123]</GONumber>
#' #    <MOC>[EXP]</MOC>
#' #  </Protein>
#' #  <Protein>
#' #    <ProteinName>4</ProteinName>
#' #    <ProteinNumber>4</ProteinNumber>
#' #    <GONumber>[GO:123123123123]</GONumber>
#' #    <MOC>[EXP]</MOC>
#' #  </Protein>
#' #  <Protein>
#' #    <ProteinName>5</ProteinName>
#' #    <ProteinNumber>5</ProteinNumber>
#' #    <GONumber>[GO:123123123123]</GONumber>
#' #    <MOC>[EXP]</MOC>
#' #  </Protein>
#' #</Family>
write_pli <- function(
  family_id, protein_name, protein_number, go_number, moc = "EXP",
  file = "") {
  
  dat <- split(data.frame(
    pname = as.character(protein_name),
    pnum  = as.character(protein_number),
    gonum = as.character(go_number),
    moc_  = as.character(moc)
  ), protein_name)
  
  dat <- sapply(dat, function(d) {
    ans <- with(d, {
      c(
        sprintf("    <ProteinName>%s</ProteinName>", unique(pname)),
        sprintf("    <ProteinNumber>%s</ProteinNumber>", unique(pnum)),
        sprintf("    <GONumber>[%s]</GONumber>", paste(gonum, collapse=", ")),
        sprintf("    <MOC>[%s]</MOC>", paste(moc_, collapse=", "))
      )
    })
    
    sprintf("  <Protein>\n%s\n  </Protein>", paste(ans, collapse = "\n"))
    
  })
  
  dat <- paste(dat, collapse = "\n")
  
  cat(paste(
    sprintf('<?xml version="1.0"?>\n<Family>\n  <FamilyID>%s</FamilyID>', family_id),
    dat, "</Family>", sep = "\n"
  ), file = file, sep = "\n")
  invisible(dat)
}

