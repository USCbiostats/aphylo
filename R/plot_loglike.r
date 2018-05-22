#' Plot LogLikelihood function of experimental data
#' @param x An object of class [aphylo()]
#' @param psi_range Numeric vector of length 2. Domain of \eqn{psi}.
#' @param mu_range Numeric vector of length 2. Domain of \eqn{mu}.
#' @param Pi_range Numeric vector of length 2. Domain of \eqn{pi}.
#' @param nlevels Integer scalar. Number of levels of each parameter to create.
#' @param plotfun Function. Either [graphics::contour()],
#' @param theta Passed to `persp`.
#' @param shade Passed to `persp`.
#' @param border Passed to `persp`.
#' @param phi Passed to `persp`.
#' @param scale Passed to `persp`.
#' [graphics::persp()], or other similar function that takes at
#' least 3 parameters, `x,y,z`.
#' @param par.args List of arguments to be passed to [par()] before
#' `plotfun` is called.
#' @param ... Aditional parameters to be passed to `plotfun`.
#' @examples 
#' # Loading data
#' data(fakeexperiment)
#' data(faketree)
#' O <- new_aphylo(fakeexperiment[,2:3], faketree)
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
plot_LogLike <- function(x, ...) UseMethod("plot_LogLike")

#' @export
#' @rdname plot_LogLike
plot_LogLike.aphylo_estimates <- function(x, ...) {
  # provided <- names(list(...))
  plot_LogLike.default(x$dat, psi = x$par[1:2], mu = x$par[3:4], eta = x$par[5:6], Pi = x$par[5], ...)
}

#' @rdname plot_LogLike
#' @template parameters
#' @templateVar psi TRUE
#' @templateVar mu TRUE
#' @templateVar Pi TRUE
#' @templateVar eta TRUE
#' @export
plot_LogLike.default <- function(
  x,
  psi_range = c(0.00001, .3),
  mu_range  = c(0.00001, .3),
  Pi_range  = c(0.00001, .3),
  psi       = rep(mean(psi_range), 2),
  mu        = rep(mean(mu_range), 2),
  eta       = c(1, 1),
  Pi        = mean(Pi_range),
  nlevels   = 30,
  plotfun   = persp,
  par.args  = list(mar=c(1, 1, 1, 1), oma=c(0,0,4,0)),
  theta     = -pi*20, 
  shade     = .7,
  border    = "steelblue",
  phi       = 30,
  scale     = TRUE,
  ...
  ) {
  
  # Adjusting values
  psi_range <- range(c(psi_range, psi))
  mu_range  <- range(c(mu_range, mu))
  Pi_range  <- range(c(Pi_range, Pi))
  
  psi_range[1] <- max(.0001, psi_range[1] - .01)
  mu_range[1] <- max(.0001, mu_range[1] - .01)
  Pi_range[1] <- max(.0001, Pi_range[1] - .01)
  
  psi_range[2] <- min(1 - .0001, psi_range[2] + .01)
  mu_range[2] <- min(1 - .0001, mu_range[2] + .01)
  Pi_range[2] <- min(1 - .0001, Pi_range[2] + .01)
  
  
  # Creating space
  PSI <- seq(psi_range[1], psi_range[2], length.out = nlevels)
  MU  <- seq(mu_range[1], mu_range[2], length.out = nlevels)
  PI  <- seq(Pi_range[1], Pi_range[2], length.out = nlevels)
  
  # Computing the actual loglike value at the selected point
  ll <- LogLike(
    x,
    psi = psi, 
    mu  = mu,
    eta = eta,
    Pi  = Pi,
    verb_ans = FALSE,
    check_dims = FALSE
  )$ll
  
  # Check PSI
  psi_z <- matrix(nrow = nlevels, ncol = nlevels)
  for (i in 1:nlevels)
    for (j in 1:nlevels)
      psi_z[i, j] <- LogLike(
        x,
        psi = c(PSI[i], PSI[j]),
        mu  = mu,
        eta = eta,
        Pi  = Pi,
        verb_ans = FALSE, 
        check_dims  = FALSE
      )$ll
  

  pi_z <- matrix(nrow = nlevels, ncol = nlevels)
  for (i in 1:nlevels)
    for (j in 1:nlevels) {
      pi_z[i, j] <- LogLike(
        x,
        psi = psi,
        mu  = c(MU[j], mu[2]),
        eta = eta,
        Pi  = PI[i],
        verb_ans = FALSE, 
        check_dims  = FALSE
      )$ll
    }
  
  mu_z <- matrix(nrow = nlevels, ncol = nlevels)
  for (i in 1:nlevels)
    for (j in 1:nlevels)
      mu_z[i, j] <- LogLike(
        x,
        psi = psi,
        mu  = c(MU[i], MU[j]),
        eta = eta,
        Pi  = Pi,
        verb_ans = FALSE, 
        check_dims  = TRUE
      )$ll
  
  # Plotting
  oldpar <- par(no.readonly = TRUE)
  if (!("mar" %in% names(par.args)))
    par.args$mar <- oldpar$mar * c(1, 1, 0.25, 0.25)
  
  # Calling plot parameters
  do.call(par, c(par.args, list(mfrow = c(2, 2))))
  
  # Actual plotting
  smu <- sprintf("%0.4f", mu)
  spsi <- sprintf("%0.4f", psi)
  sPi <- sprintf("%0.4f", Pi)
  
  # Replacing the infs
  psi_z[is.infinite(psi_z)] <- min(psi_z[!is.infinite(psi_z)])
  mu_z[is.infinite(mu_z)]   <- min(mu_z[!is.infinite(mu_z)])
  pi_z[is.infinite(pi_z)]   <- min(pi_z[!is.infinite(pi_z)])
  
  pmat <- plotfun(
    PSI,
    PSI,
    psi_z,
    xlab = expression(psi[0]),
    ylab = expression(psi[1]),
    main = bquote(mu[0] == .(smu[1])~ mu[1] == .(smu[2]) ~ and ~ pi == .(sPi)),
    theta = theta, shade = shade, border = border, phi = phi, scale = scale,
    ...
  )
  
  # Drawing the point
  if (identical(plotfun, persp))
    points(trans3d(psi[1], psi[2], ll, pmat), col = "black", pch = 25, bg = "red", cex = 1.5)
  else 
    points(psi[1], psi[2], col = "black", pch = 25, bg = "red", cex = 1.5)
  
  
  pmat <- plotfun(
    MU,
    MU,
    mu_z,
    xlab = expression(mu[0]),
    ylab = expression(mu[1]),
    main = bquote(psi[0] == .(spsi[1]) ~ psi[1] == .(spsi[2]) ~ and ~ pi == .(sPi)),
    theta = theta, shade = shade, border = border, phi = phi, scale = scale,
    ...
  )
  
  # Drawing the point
  if (identical(plotfun, persp))
    points(trans3d(mu[1], mu[2], ll, pmat), col = "black", pch = 25, bg = "red", cex = 1.5)
  else
    points(mu[1], mu[2], col = "black", pch = 25, bg = "red", cex = 1.5)
  
  pmat <- plotfun(
    PI,
    MU,
    pi_z,
    xlab = expression(pi),
    ylab = expression(mu[0]),
    main = bquote(psi[0] == .(spsi[1])~ psi[1]==.(spsi[2]) ~ and ~ mu[1] == .(smu[2])),
    theta = theta, shade = shade, border = border, phi = phi, scale = scale,
    ...
  )
  
  # Drawing the point
  if (identical(plotfun, persp))
    points(trans3d(Pi, mu[1], ll, pmat), col = "black", pch = 25, bg = "red", cex = 1.5)
  else
    points(Pi, mu[1], col = "black", pch = 25, bg = "red", cex = 1.5)
  
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
      Pi  = list(PI, pi_z),
      eta = eta
    )
  )
}


#' Multiavariate plot (surface)
#' @param fun A function that receives 2 or more parameters and returns a single
#' number.
#' @param params Numeric vector with the default parameters.
#' @param domain (optional) Named list with as many elements as parameters. Specifies the
#' domain of the function.
#' @param nlevels Integer. Number of levels.
#' @param args List of named arguments to be passed to `fun`.
#' @param plotfun Function that will be used to plot `x,y,z`.
#' @param plot Logical. When `FALSE` skips plotting.
#' @param postplot Function to be called after `plotfun`. Should recieve a vector
#' with the current parameters.
#' @param ... Further arguments passed to `plotfun`.
#' @export
#' 
#' @examples 
#' # Example: A model with less parameters
#' x <- sim_annotated_tree(20)
#' ans <- aphylo_mcmc(x ~ psi, control=list(nbatch=1e4))
#' 
#' # Creating the multivariate plot (using by default image)
#' plot_multivariate(
#'   function(...) {
#'     ans$fun(unlist(list(...)), priors = ans$priors, dat = ans$dat, verb_ans = FALSE)
#'   },
#'   params = ans$par
#' )
plot_multivariate <- function(
  fun,
  params,
  domain,
  nlevels = 20,
  args    = list(),
  plotfun = image,
  plot    = TRUE,
  postplot = function(params) {
    points(params, cex = 2, pch=3, col="red")
  },
  ...
) {
  
  # How many params
  k <- length(params)
  
  # Specifying names
  if (!length(names(params)))
    names(params) <- sprintf("par%02i", 1L:k)
  pnames <- names(params)
  
  # Specifying limits
  if (missing(domain)) 
    domain <- structure(
      lapply(pnames, function(i) c(1e-10, 1 - 1e-10)),
      names = pnames
      )
  
  # Listing how many combs
  plots <- utils::combn(pnames, 2)
  
  ans <- vector("list", ncol(plots))
  for (p in seq_len(ncol(plots))) {
    
    # Computing points
    z <- matrix(ncol=nlevels, nrow=nlevels)
    
    x <- seq(domain[[plots[1,p]]][1], domain[[plots[1, p]]][2], length.out = nlevels)
    y <- seq(domain[[plots[2,p]]][1], domain[[plots[2, p]]][2], length.out = nlevels)
    
    tmppar <- as.list(params)
    for (i in 1:nlevels)
      for (j in 1:nlevels) {
        
        # Replacing the values
        tmppar[[plots[1, p]]] <- x[i]
        tmppar[[plots[2, p]]] <- y[j]
        
        # Computing the z value
        z[i, j] <- do.call(fun, c(tmppar, args))
        
      }
       
    ans[[p]] <- list(
      x = x,
      y = y,
      z = z,
      xlab = plots[1, p],
      ylab = plots[2, p]
      )
    
  }
  
  if (plot) {
    op <- graphics::par(mfrow = grDevices::n2mfrow(ncol(plots)))
    on.exit(graphics::par(op))
    for (p in seq_along(ans)) {
      
      do.call(plotfun, c(ans[[p]], list(...)))
      postplot(params[plots[,p,drop=TRUE]])
      
      
    }
  }
  
    
  invisible(ans)
  
}



