#' Accuracy calculation as defined in Engelhardt et al. (2011)
#' 
#' Uses SIFTER's 2011 definition of accuracy, where a protein is tagged as
#' accurately predicted if the highest ranked prediction matches it.
#' 
#' @param pred A matrix of predictions, or an [aphylo_estimates] object.
#' @param lab A matrix of labels (0,1,NA, or 9 if `nine_na = TRUE`).
#' @param tol Numeric scalar. Predictions within `tol` of the max score
#' will be tagged as the prediction made by the model (see deails).
#' @param highlight Pattern passed to [sprintf] used to highlight
#' predicted functions that match the observed.
#' @param nine_na Treat 9 as NA.
#' @param ... Further arguments passed to the method. In the case of `aphylo_estimates`,
#' the arguments are passed to [predict.aphylo_estimates()].
#' 
#' @details 
#' The analysis is done at the protein level. For each protein, the function
#' compares the YES annotations of that proteins with the predicted by the model.
#' The algorithm selects the predicted annotations as those that are within
#' `tol` of the maximum score.
#' 
#' This algorithm doesn't take into account NOT annotations (0s), which are
#' excluded from the analysis.
#' 
#' When `highlight = ""`, no highlight is done.
#' @returns 
#' A data frame with `Ntip()` rows and four variables. The variables are:
#' - Gene: Label of the gene
#' - Predicted: The assigned gene function.
#' - Observed: The true set of gene functions.
#' - Accuracy: The measurement of accuracy according to Engelhardt et al. (2011).
#' 
#' @examples 
#' set.seed(81231)
#' atree <- raphylo(50, psi = c(0,0), P = 3)
#' ans <- aphylo_mcmc(atree ~ mu_d + mu_s + Pi)
#' 
#' accuracy_sifter(ans)
#' 
#' @export
accuracy_sifter <- function(
  pred,
  lab,
  tol       = 1e-10,
  highlight = "",
  ...) UseMethod("accuracy_sifter")

#' @export
#' @rdname accuracy_sifter
accuracy_sifter.aphylo_estimates <- function(
  pred, lab, tol = 1e-10,
  highlight = "", ...
  ) {
  
  if (Ntrees(pred) > 1)
    stop("No support for multiple trees (yet).")
  
  ids <- 1:Ntip(pred)
  
  accuracy_sifter(
    pred    = predict(pred, ids = list(ids),...)[ids,,drop=FALSE],
    lab     = pred$dat$tip.annotation,
    tol     = tol,
    nine_na = TRUE,
    highlight = ""
  )
  
}

#' @export
#' @rdname accuracy_sifter
accuracy_sifter.default <- function(
  pred, lab, tol = 1e-10,
  highlight = "",
  nine_na = TRUE, ...
  ) {
   
  # Coercing to the right type
  pred <- as.matrix(pred)
  lab  <- as.matrix(lab)
  
  if (!all(dim(pred) == dim(lab)))
    stop("All the dimensions should match.")
  
  # Identifying observations
  if (nine_na) {
    
    # All in the right place?
    if (!all(as.vector(lab) %in% c(0,1,9)))
        stop("The matrix -lab- has values other than 0, 1, or 9.")
        
    # Looking for values that are not missing
    lab[lab == 0] <- 9
    ids <- which(rowSums(lab) != 9*ncol(lab))
    
  } else {
    
    # All in the right place?
    if (!all(as.vector(lab) %in% c(0,1,NA)))
        stop("The matrix -lab- has values other than 0, 1, or NA.")
        
    lab[is.na(lab) | lab == 0] <- 9
    # Looking for values that are not missing
    ids <- which(rowSums(lab) != 9*ncol(lab))
    
  }
  
  # Getting the function names
  fnames <- colnames(pred)
  
  ans <- lapply(ids, function(i) {
    top_pred <- fnames[which(abs(pred[i,] - max(pred[i,])) < tol)]
    true_ann <- fnames[which(lab[i,] == 1)]
    
    its_a_match <- which(top_pred %in% true_ann)
    if (length(its_a_match) && highlight != "")
      top_pred[its_a_match] <- sprintf(highlight, top_pred[its_a_match])
    
    data.frame(
      Predicted = paste(top_pred, collapse=","),
      Observed  = paste(true_ann, collapse=","),
      Accuracy  = ifelse(length(its_a_match), 1, 0)
    )  
  })
  
  cbind(Gene = rownames(pred)[ids],do.call(rbind, ans))
}

