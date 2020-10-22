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
#' @export
write_pli <- function(family_id, protein_name, protein_number, go_number, moc = "EXP") {
  
  dat <- split(cbind(protein_name, go_number, moc), protein_name)
  
  dat <- sapply(dat, function(d) {
    ans <- with(d, {
      c(
        sprintf("\t\t<ProteinName>%s</ProteinName>", unique(protein_name)),
        sprintf("\t\t<ProteinNumber>%s</ProteinNumber>", unique(protein_number)),
        sprintf("\t\t<GONumber>[%s]</GONumber>", paste(go_number, collapse=", ")),
        sprintf("\t\t<MOC>[%s]</MOC>", paste(moc, collapse=", ")),
      )
    })
    
    sprintf("\t<Protein>\n%s\n\t</Protein>", paste(ans, collapse = "\n"))
    
  })
  
  dat <- paste(dat, collapse = "\n")
  
  paste(
    sprintf('<?xml version="1.0"?><Family>\n\t<FamilyID>%s</FamilyID>', family_id),
    dat, "</Family>", sep = "\n"
  )
  
}

