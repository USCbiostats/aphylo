#' Plot Log-Likelihood function of the model
#' @param x An object of class [aphylo()]
#' @inheritParams plot_multivariate
#' @param ... Aditional parameters to be passed to `plotfun`.
#' @examples 
#' # Loading data
#' data(fakeexperiment)
#' data(faketree)
#' O <- new_aphylo(fakeexperiment[,2:3], faketree)
#' 
#' # Baseline plot (all parameters but Pi)
#' plot_logLik(O)
#' 
#' # No psi parameter
#' plot_logLik(O ~ mu + Pi + eta)
#' 
#' @export

#' @export
#' @rdname plot_logLik
plot_logLik <- function(x, sets, ...) UseMethod("plot_logLik")

#' @export
#' @rdname plot_logLik
plot_logLik.aphylo <- function(x, sets, ...) {
  
  # Creating formula
  env <- parent.frame()
  if (inherits(x, "formula"))
    model <- aphylo_formula(x, env = env)
  else {
    x <- substitute(x ~ psi + mu + eta, list(x = substitute(x)))
    x <- stats::as.formula(x, env)
    model <- aphylo_formula(x, env = env)
  }
  
  plot_logLik.aphylo_estimates(
    list(
      dat    = model$dat,
      params = model$params,
      priors = function(p) 1,
      fun    = model$fun,
      ll     = model$fun(model$params, function(p) 1, model$dat, FALSE)
    ),
    sets = sets,
    ...
  )
  
}

#' @export
#' @rdname plot_logLik
plot_logLik.formula <- plot_logLik.aphylo

#' @export
#' @rdname plot_logLik
plot_logLik.aphylo_estimates <- function(x, sets,...) {
  
  # Collecting dots
  dots <- list(...)
  
  # if (!length(dots))
  
  # Generating combinations
  if (missing(sets)) {
    
    sets <- NULL
    for (p in c("psi", "mu", "eta"))
      if (any(grepl(p, names(x$par))))
        sets <- cbind(sets, paste0(p, 0:1))

  }
  
  # Finding arrangement
  mfrow <- if(ncol(sets) == 1)
    c(1,1)
  else if (ncol(sets) == 2)
    c(1, 2)
  else if (ncol(sets) >= 3)
    c(2, 2)
  
  # Plotting
  op <- graphics::par(mar = c(1, 1, .2, .2), oma = c(0, 0, 2, 0))
  on.exit(par(op))
  plot_multivariate(
    function(...) {
      x$fun(unlist(list(...)), priors = x$priors, dat = x$dat, verb_ans = FALSE)
    },
    params  = x$par,
    sets    = sets,
    plotfun = function(...) {
      
      # Capturing arguments
      dots <- list(...)
      nrz <- nrow(dots$z)
      ncz <- ncol(dots$z)
      
      # Fixing ranges
      dots$zlim     <- as.vector(stats::quantile(dots$z, c(.025, 1)))
      dots$ticktype <- "detailed"
      dots$z[dots$z > dots$zlim[2]] <- dots$zlim[2]
      dots$z[dots$z < dots$zlim[1]] <- dots$zlim[1]
      
      # Creating colors
      nbcol       <- 100
      color       <- viridis::viridis(nbcol)
      zfacet      <- dots$z[-1, -1] + dots$z[-1, -ncz] + dots$z[-nrz, -1] + dots$z[-nrz, -ncz]
      facetcol    <- cut(zfacet, nbcol)
      dots$col    <- color[facetcol]
      dots$border <- grDevices::adjustcolor(dots$col, red.f = .5, green.f = .5, blue.f = .5)
      

      # Making some room for the labels
      dots$xlab <- paste0("\n\n", dots$xlab)
      dots$ylab <- paste0("\n\n", dots$ylab)
      
      # Creating colors
      do.call(graphics::persp, dots)
    },
    theta   = -pi*20, 
    shade   = .7,
    zlab    = "",
    mfrow   = mfrow,
    phi     = 30,
    postplot = function(par, res) {
      graphics::points(
        grDevices::trans3d(par[1], par[2], x$ll, res),
        col = "black", pch = 25, bg = "red", cex = 1.5
        )
    },
    ...
  )
  
  graphics::par(mfrow = c(1,1), xpd=NA)
  graphics::mtext(paste0("Log L(", paste0(names(x$par), collapse = ","), ")"),
        outer=TRUE)
  
  # For later on
  # pmat <- plotfun(
  #   MU,
  #   MU,
  #   mu_z,
  #   xlab = expression(mu[0]),
  #   ylab = expression(mu[1]),
  #   main = bquote(psi[0] == .(spsi[1]) ~ psi[1] == .(spsi[2]) ~ and ~ pi == .(sPi)),
  #   theta = theta, shade = shade, border = border, phi = phi, scale = scale,
  #   ...
  # )
  
  
  
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
#' @param sets (optional) Character matrix of size `2 x # of combinations`.
#' contains the names of the pairs to plot. If nothing passed, the function will
#' generate all possible combinations as `combn(names(params), 2)`.
#' @param mfrow Passed to [graphics::par].
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
#'   sets = matrix(c("mu0", "mu1", "psi0", "psi1"), ncol=2),
#'   params = ans$par
#' )
plot_multivariate <- function(
  fun,
  params,
  domain,
  sets,
  nlevels = 20,
  args    = list(),
  plotfun = graphics::image,
  plot    = TRUE,
  postplot = function(params, res) {
    points(params, cex = 2, pch=3, col="red")
  },
  mfrow   = NULL,
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
  if (missing(sets))
    sets <- utils::combn(pnames, 2)
  
  ans <- vector("list", ncol(sets))
  for (p in seq_len(ncol(sets))) {
    
    # Computing points
    z <- matrix(ncol=nlevels, nrow=nlevels)
    
    x <- seq(domain[[sets[1,p]]][1], domain[[sets[1, p]]][2], length.out = nlevels)
    y <- seq(domain[[sets[2,p]]][1], domain[[sets[2, p]]][2], length.out = nlevels)
    
    tmppar <- as.list(params)
    for (i in 1:nlevels)
      for (j in 1:nlevels) {
        
        # Replacing the values
        tmppar[[sets[1, p]]] <- x[i]
        tmppar[[sets[2, p]]] <- y[j]
        
        # Computing the z value
        z[i, j] <- do.call(fun, c(tmppar, args))
        
      }
    
    # Storing the results
    ans[[p]] <- list(
      x = x,
      y = y,
      z = z,
      xlab = sets[1, p],
      ylab = sets[2, p]
      )
    
  }
  
  # Should we plot?
  if (plot) {
    
    if (!length(mfrow))
      mfrow <- grDevices::n2mfrow(ncol(sets))
    
    op <- graphics::par(mfrow = mfrow, xpd=NA)
    on.exit(graphics::par(op))
    for (p in seq_along(ans)) {
      
      res <- do.call(plotfun, c(ans[[p]], list(...)))
      postplot(params[sets[,p,drop=TRUE]], res)
      
      
    }
  }
  
    
  invisible(ans)
  
}



