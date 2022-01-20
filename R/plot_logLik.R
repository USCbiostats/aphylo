#' Plot Log-Likelihood function of the model
#' @param x An object of class [aphylo()]
#' @inheritParams plot_multivariate
#' @param ... Aditional parameters to be passed to `plotfun`.
#' @returns 
#' NULL (invisible). Generates a plot of the loglikelihood of the model.
#' @examples 
#' # Loading data
#' data(fakeexperiment)
#' data(faketree)
#' O <- new_aphylo(fakeexperiment[,2:3], tree = as.phylo(faketree))
#' 
#' # Baseline plot (all parameters but Pi)
#' plot_logLik(O)
#' 
#' # No psi parameter
#' plot_logLik(O ~ mu_d + Pi + eta)
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
    x <- substitute(x ~ psi + mu_d + mu_s + eta, list(x = substitute(x)))
    x <- stats::as.formula(x, env)
    model <- aphylo_formula(x, env = env)
  }
  
  plot_logLik.aphylo_estimates(
    list(
      dat    = model$dat,
      params = model$params,
      priors = function(p) 1,
      fun    = model$fun,
      ll     = model$fun(model$params, model$dat, function(p) 1, FALSE)
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
    for (p in c("psi", "mu_d", "eta"))
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
  dat0 <- new_aphylo_pruner(x$dat)
  plot_multivariate(
    function(...) {
      x$fun(unlist(list(...)), priors = x$priors, dat = dat0, verb_ans = FALSE)
    },
    params  = x$par,
    sets    = sets,
    plotfun = function(...) {
      
      # Capturing arguments
      dots <- list(...)
      nrz <- nrow(dots$z)
      ncz <- ncol(dots$z)
      
      # Fixing ranges
      dots$zlim     <- as.vector(stats::quantile(dots$z, c(.025, 1), na.rm=TRUE)) + c(-.01, +.01)
      dots$ticktype <- "detailed"
      dots$z[dots$z > dots$zlim[2]] <- dots$zlim[2]
      dots$z[dots$z < dots$zlim[1]] <- dots$zlim[1]
      
      # Creating colors
      nbcol       <- 100
      color       <- # viridisLite::viridis(100)
        c("#440154FF", "#450558FF", "#46085CFF", "#470D60FF", "#471063FF", 
          "#481467FF", "#481769FF", "#481B6DFF", "#481E70FF", "#482173FF", 
          "#482576FF", "#482878FF", "#472C7AFF", "#472F7CFF", "#46327EFF", 
          "#453581FF", "#453882FF", "#443B84FF", "#433E85FF", "#424186FF", 
          "#404587FF", "#3F4788FF", "#3E4A89FF", "#3D4D8AFF", "#3C508BFF", 
          "#3B528BFF", "#39558CFF", "#38598CFF", "#375B8DFF", "#355E8DFF", 
          "#34608DFF", "#33638DFF", "#32658EFF", "#31688EFF", "#2F6B8EFF", 
          "#2E6D8EFF", "#2D708EFF", "#2C718EFF", "#2B748EFF", "#2A768EFF", 
          "#29798EFF", "#287C8EFF", "#277E8EFF", "#26818EFF", "#26828EFF", 
          "#25858EFF", "#24878EFF", "#238A8DFF", "#228D8DFF", "#218F8DFF", 
          "#20928CFF", "#20938CFF", "#1F968BFF", "#1F998AFF", "#1E9B8AFF", 
          "#1F9E89FF", "#1FA088FF", "#1FA287FF", "#20A486FF", "#22A785FF", 
          "#24AA83FF", "#25AC82FF", "#28AE80FF", "#2BB07FFF", "#2EB37CFF", 
          "#31B67BFF", "#35B779FF", "#39BA76FF", "#3DBC74FF", "#41BE71FF", 
          "#47C06FFF", "#4CC26CFF", "#51C56AFF", "#56C667FF", "#5BC863FF", 
          "#61CA60FF", "#67CC5CFF", "#6DCD59FF", "#73D056FF", "#78D152FF", 
          "#7FD34EFF", "#85D54AFF", "#8CD646FF", "#92D741FF", "#99D83DFF", 
          "#A0DA39FF", "#A7DB35FF", "#ADDC30FF", "#B4DE2CFF", "#BBDE28FF", 
          "#C2DF23FF", "#C9E020FF", "#D0E11CFF", "#D7E219FF", "#DDE318FF", 
          "#E4E419FF", "#EBE51AFF", "#F1E51DFF", "#F7E620FF", "#FDE725FF"
        )
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
  #   mu_d,
  #   mu_d,
  #   mu_z,
  #   xlab = expression(mu_d[0]),
  #   ylab = expression(mu_d[1]),
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
#' @returns 
#' A list of length `length(sets)`, each with the following:
#' - x,y,z vectors of coordinates. 
#' - xlab,ylab vectors with the corresponding labels.
#' 
#' @examples 
#' # Example: A model with less parameters
#' set.seed(1231)
#' x <- raphylo(20)
#' ans <- aphylo_mcmc(
#'   x ~ psi + mu_d + mu_s,
#'   control = list(nsteps = 1e3, burnin = 0)
#'   )
#' 
#' # Creating the multivariate plot (using by default image)
#' plot_multivariate(
#'   function(...) {
#'     ans$fun(unlist(list(...)), priors = ans$priors, dat = ans$dat, verb_ans = FALSE)
#'   },
#'   sets = matrix(c("mu_d0", "mu_d1", "psi0", "psi1"), ncol=2),
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
      lapply(pnames, function(i) c(1e-5, 1 - 1e-5)),
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
    
    # We don't want values that are too big
    z[abs(z) > 1e100] <- NA
    
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



