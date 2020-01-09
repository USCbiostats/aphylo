library(aphylo)

#' This is a very simple example, here we only condition on:
#' - # changes
#' - max(blength)
#' - either or (0,0)->(0,1) (0,0)->(1,0), (1,1) -> (1,0) (1,1) -> (0,1)
pr <- function(x_off, x_par, mu, blength = NULL) {
  
  # Scalars
  P <- length(x_off)
  # if (is.null(blength))
    blength <- rep(1, P)
  
  # Computing powerset offspring
  pset <- aphylo:::states(P)
  
  # Calculating suff stats
  
  s_nfuns   <- colSums(t(pset) != x_par)
  s_xmax    <- apply(t(pset)*blength, 2, max)
  s_eor     <- as.integer(colSums(abs(t(pset) - x_par)) == 1L)
  
  s_nfuns_i <- sum(x_off != x_par)
  s_xmax_i  <- max(x_off * blength)
  s_eor_i   <- as.integer(sum(abs(x_off - x_par)) == 1)
  
  S  <- sum(exp(cbind(s_nfuns, s_xmax, s_eor) %*% mu))
  Si <- exp(sum(c(s_nfuns_i, s_xmax_i, s_eor_i) * mu))
  
  list(
    p = Si/S,
    stats = cbind(s_nfuns, s_xmax, s_eor)
  )
  
}

# All cases printer
printer <- function(x0, x1, p, mu, stats) {

  cat(
    "\\pagebreak\n\n\\noindent {\\bf \\large Case where $X_{par} = (", 
    paste0(sprintf("%d", x0), collapse=", ")
    ,")$}\n\n", sep="")
  
  cat(
    "In this case the parameters are $\\mu = (",
    paste0(sprintf("%.2f", mu), collapse = ", "),
    ")$. We are estimating ", length(x0), " genes jointly for which SIFTER ",
    "would need ", 2^(2*length(x0)), " parameters.",
    "\n\n",
    sep=""
    )
  
  cat("\\begin{table}[!h]\n\\begin{tabular}{lcr}\n\\toprule\n")  
  
  cat("Transition & $\\mbox{(n changes, max blength, changes = 1)}$ & Prob\\\\\n\\midrule")
  
  if (length(x0) == 2) {
    cat(sprintf(
      "$(%d, %d) \\to (%d, %d)$ & (%.2f, %.2f) &  %.2f \\\\\n",
      x0[1], x0[2],
      x1[,1], x1[2],
      stats[,1], stats[,3],
      p
      ), sep = "")
  } else if (length(x0) == 3) {
    cat(sprintf(
      "$(%d, %d, %d) \\to (%d, %d, %d)$ & (%.2f,%.2f) &  %.2f \\\\\n",
      x0[1], x0[2], x0[3],
      x1[,1], x1[,2], x1[,3],
      stats[,1], stats[,3],
      p
    ), sep = "")
  }
  cat("\\bottomrule\\end{tabular}\n\\end{table}\n\n")
  
  cat("The marginal probabilities in this case are:\n $(",
      paste0(sprintf("%.2f", colSums(x1*p)), collapse=", "),
      ")$.\n\n", sep=""
      )
}


# Case 1: Moving from c(1,1,1) -------------------------------------------------

L <- list(c(1,1,1), c(1,1,0), c(1,0,0), c(0,0,0))

for (x0 in L) {

  bl <- c(1, 2, 1.5)
  mu <- c(-1, 1, 1) 
  
  # Cases
  allcases <- aphylo:::states(length(x0))
  ans <- vector("list",nrow(allcases))
  for (i in seq_along(ans))
    ans[[i]] <- pr(x_off = allcases[i,], x_par = x0, mu = mu, blength = bl)

  # print(sum(ans))
  
  printer(x0, allcases, sapply(ans, "[[", "p"), mu, ans[[1]]$stats)
}