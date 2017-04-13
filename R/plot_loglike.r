#' Plot LogLikelihood function of experimental data
#' @param x An object of class \code{\link[phylogenetic:new_aphylo]{aphylo}}
#' @param psi_range Numeric vector of length 2. Domain of \eqn{psi}.
#' @param mu_range Numeric vector of length 2. Domain of \eqn{mu}.
#' @param Pi_range Numeric vector of length 2. Domain of \eqn{pi}.
#' @param nlevels Integer scalar. Number of levels of each parameter to create.
#' @param plotfun Function. Either \code{\link[graphics:contour]{contour}},
#' \code{\link[graphics:persp]{persp}}, or other similar function that takes at
#' least 3 parameters, \code{x,y,z}.
#' @param par.args List of arguments to be passed to \code{\link{par}} before
#' \code{plotfun} is called.
#' @param ... Aditional parameters to be passed to \code{plotfun}.
#' @examples 
#' # Loading data
#' data(fakeexperiment)
#' data(faketree)
#' O <- new_aphylo(fakeexperiment, faketree[,c("ParentId", "NodeId")], "LeafId")
#' 
#' # Nice personalized plot
#' plot_LogLike(O, nlevels = 60, plotfun = persp, theta = -pi*20, 
#'   shade=.7, border="darkblue", phi=30, scale=TRUE, 
#'   par.args = list(mar=c(1, 1, 1, 1), oma=c(0,0,4,0)),
#'   ticktype = "detailed")
#'   
#' # Adding title
#' mtext(
#'   "LogLikelihood",
#'   side=3, outer=FALSE, line = 2, cex=1.25)
#' @export
plot_LogLike <- function(
  x,
  psi_range = c(0.00001, .1),
  mu_range  = c(0.00001, .9999),
  Pi_range  = c(0.00001, .9999),
  nlevels   = 30,
  plotfun   = persp,
  par.args  = list(),
  ...
  ) {
  
  # Creating space
  PSI <- seq(psi_range[1], psi_range[2], length.out = nlevels)
  MU  <- seq(mu_range[1], mu_range[2], length.out = nlevels)
  PI  <- seq(Pi_range[1], Pi_range[2], length.out = nlevels)
  
  # Check PSI
  mu    <- mean(MU)
  Pi    <- mean(PI)
  psi_z <- matrix(nrow = nlevels, ncol = nlevels)
  for (i in 1:nlevels)
    for (j in 1:nlevels)
      psi_z[i, j] <- LogLike(
        x$annotations,
        x$offspring,
        x$noffspring,
        c(PSI[i], PSI[j]),
        c(mu, mu),
        c(Pi, 1 - Pi),
        verb_ans = FALSE
      )$ll
  
  psi  <- mean(PSI)
  pi_z <-
    matrix(nrow = nlevels, ncol = nlevels) # vector("numeric", nlevels)
  for (i in 1:nlevels)
    for (j in 1:nlevels) {
      pi_z[i, j] <- LogLike(
        x$annotations,
        x$offspring,
        x$noffspring,
        c(psi, psi),
        c(MU[j], mu),
        c(PI[i], 1 - PI[i]),
        verb_ans = FALSE
      )$ll
    }
  
  mu_z <- matrix(nrow = nlevels, ncol = nlevels)
  for (i in 1:nlevels)
    for (j in 1:nlevels)
      mu_z[i, j] <- LogLike(
        x$annotations,
        x$offspring,
        x$noffspring,
        c(psi, psi),
        c(MU[i], MU[j]),
        c(Pi, 1 - Pi),
        verb_ans = FALSE
      )$ll
  
  # Plotting
  oldpar <- par(no.readonly = TRUE)
  if (!("mar" %in% names(par.args)))
    par.args$mar <- oldpar$mar * c(1, 1, 0.25, 0.25)
  
  # Calling plot parameters
  do.call(par, c(par.args, list(mfrow = c(2, 2))))
  
  # Actual plotting
  mu <- sprintf("%0.4f", mu)
  psi <- sprintf("%0.4f", psi)
  Pi <- sprintf("%0.4f", Pi)
  
  plotfun(
    PSI,
    PSI,
    psi_z,
    xlab = expression(psi[0]),
    ylab = expression(psi[1]),
    main = bquote(mu == .(mu) ~ and ~ pi == .(Pi)),
    ...
  )
  
  plotfun(
    MU,
    MU,
    mu_z,
    xlab = expression(mu[0]),
    ylab = expression(mu[1]),
    main = bquote(psi == .(psi) ~ and ~ pi == .(Pi)),
    ...
  )
  
  plotfun(
    PI,
    MU,
    pi_z,
    xlab = expression(pi),
    ylab = expression(mu[0]),
    main = bquote(psi == .(psi) ~ and ~ mu[1] == .(mu)),
    ...
  )
  
  # Adding legend
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  legend(
    "center",
    legend = expression(
      pi ~ Root ~ node ~ probabilities,
      psi ~ Misclassification ~ probabilities,
      mu ~ Loss / Gain ~ probabilities
    ),
    bty = "n"
  )
  
  # Restoring parameters
  par(oldpar)
  
  invisible(
    list(
      psi = list(PSI, psi_z),
      mu  = list(MU, mu_z),
      Pi  = list(PI, pi_z)
    )
  )
}

