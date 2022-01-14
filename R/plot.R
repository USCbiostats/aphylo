#' Set of colors from dput(RColorBrewer::brewer.pal(7, "RdBu"), file = "")
#' @noRd
.aphyloColors <- c("#B2182B", "#EF8A62", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#67A9CF", 
                   "#2166AC")


#' Plot and print methods for `aphylo` objects
#' 
#' @param x An object of class `aphylo`.
#' @param y Ignored.
#' @param ... Further arguments passed to [ape::plot.phylo].
#' @param prop Numeric scalar between 0 and 1. Proportion of the device that the
#' annotations use in `plot.aphylo`.
#' @param rect.args List of arguments passed to [graphics::rect].
#' @name aphylo-methods
#' @param node.type.col,node.type.size Vectors of length 2. In the case of 
#' `node.type.col` the color of the duplication and other nodes. `node.type.size`
#' sets the size of circles.
#' @param as_ci Integer vector. Internal use only.
#' @details The `plot.aphylo` function is a wrapper of [ape::plot.phylo].
#' 
#' @export
#' @return In the case of `plot.aphylo`, `NULL`.
#' @family aphylo methods
#' @export
#' @examples 
#' set.seed(7172)
#' atree <- raphylo(20)
#' plot(atree)
plot.aphylo <- function(
  x,
  y              = NULL,
  prop           = .15, 
  node.type.col  = c(dupl = "black", other = "gray"),
  node.type.size = c(dupl = 0, other = 0), 
  rect.args      = list(),
  as_ci          = NULL,
  ...
) {
  
  # Coercing into phylo
  phylo <- as.phylo(x)
  dots  <- list(...)
  
  # Some defaults
  if (!length(dots$cex))
    dots$cex <- .75
  if (!length(dots$show.node.label))
    dots$show.node.label <- FALSE
  if (!length(dots$font))
    dots$font <- 1
  if (!length(dots$main))
    dots$main <- "Annotated Phylogenetic Tree"
  if (!length(dots$align.tip.label))
    dots$align.tip.label <- TRUE
  
  if (length(dots$type) && dots$type != "phylogram")
    stop("Only `phylograph` is currently supported.")
  
  # Size of the device
  dev_size <- graphics::par("din")
  
  # How much space for the annotations
  labwidth <- dev_size[1]*prop
  
  op <- graphics::par(
    mai = graphics::par("mai")*c(0, 1, 1, 0) + c(labwidth, 0, 0, labwidth)
  )
  
  on.exit(graphics::par(op))
  do.call(graphics::plot, c(list(x=phylo), dots))
  
  # Capturing the parameters from the `ape` package
  plot_pars <- utils::getFromNamespace(".PlotPhyloEnv", "ape")
  
  # Adding node type
  nodes <- (Ntip(x) + 1):Nnode(x, internal.only=FALSE)
  nodes <- with(plot_pars$last_plot.phylo, cbind(xx[nodes], yy[nodes]))
  
  tips <- with(plot_pars$last_plot.phylo, cbind(xx, yy))
  tips <- tips[1L:ape::Ntip(phylo),,drop=FALSE]
  
  # Printing node types --------------------------------------------------------  
  if (all(node.type.size == 0))
    node.type.size <- rep(diff(range(c(tips[,1], nodes[,1])))/90, 2)
  
  graphics::symbols(
    x       = nodes[,1] - node.type.size[x$node.type + 1L]*1.3,
    y       = nodes[,2],
    circles = node.type.size[x$node.type + 1L],
    add     = TRUE,
    inches  = FALSE,
    fg      = "lightgray",
    bg      = node.type.col[x$node.type + 1L]
  )
  
  # Plotting annotations (making room first) -----------------------------------
  yspacing <- range(tips[,2])
  yspacing <- (yspacing[2] - yspacing[1])/(nrow(tips) - 1)/2
  
  op2 <- graphics::par(
    mai = graphics::par("mai")*c(1,0,1,0) +
      c(0, dev_size[1]*(1-prop), 0, dev_size[1]*.025)
  )
  on.exit(graphics::par(op2), add=TRUE)
  
  graphics::plot.window(c(0, 1), range(tips[,2]), new = FALSE, xaxs = "i")

  nfun      <- ncol(x$tip.annotation)
  yran      <- range(tips[,2])
  
  # If we are plotting the CI, then it is only two blocks
  nblocks <- ifelse(length(as_ci), 2, nfun)
  ndraws  <- 1L
  for (f in 1:nfun) {
    
    rect.args$xleft   <- (ndraws - 1)/nblocks + 1/nblocks*.05
    rect.args$ybottom <- tips[,2] - yspacing
    rect.args$xright  <- ndraws/nblocks - 1/nblocks*.05
    rect.args$ytop    <- tips[,2] + yspacing
    rect.args$xpd     <- NA
    rect.args$col     <- blue(x$tip.annotation[,f])
    # rect.args$border  <- blue(x$tip.annotation[,f])
    rect.args$col[x$tip.annotation[, f] == 9L] <- "white"
    rect.args$border  <- rect.args$col
    
    if (!length(rect.args$xpd)) rect.args$xpd <- NA
    if (!length(rect.args$lwd)) rect.args$lwd<-.5
    
    # Draing a background
    graphics::rect(
      xleft   = rect.args$xleft,
      xright  = rect.args$xright,
      ytop    = max(rect.args$ytop),
      ybottom = min(rect.args$ybottom),
      col     = "darkgray",
      border  = "lightgray"
    )
    
    # In the case that the function plotted is actually a
    # confidence interval, we plot bars instead
    if (length(as_ci) && f == as_ci[1L]) {
      
      # The locations of the x coordinates change
      barwd <- rect.args$xright - rect.args$xleft
      x0    <- rect.args$xleft
      x1    <- rect.args$xright
      
      xmid  <- rect.args$xright + barwd/2
      
      rect.args$xleft  <- x0 + barwd * x$tip.annotation[, f]
      rect.args$xright <- x1 - barwd *
        (1 - x$tip.annotation[, as_ci[3L]])
      
      # Assuring a minimum
      rect.args$xleft[rect.args$xleft < (x0 + barwd*.1)] <- x0 + barwd*.1
      rect.args$xright[rect.args$xright > (x1 - barwd*.1)] <- x1 - barwd*.1
      
      # Figuring out the right color, we use the middle one
      rect.args$col    <- blue(x$tip.annotation[, as_ci[2]])
      rect.args$border <- rect.args$col
      
      rect.args$ytop    <- rect.args$ytop - yspacing*.2
      rect.args$ybottom <- rect.args$ybottom + yspacing*.2
      
    } 
    
    # Drawing rectangles
    do.call(graphics::rect, rect.args)
    
    # Adding function label
    graphics::text(
      x = (2 * ndraws - 1)/nblocks/2 - 1/nblocks/2,
      y = yran[1] - graphics::strheight(
        colnames(x$tip.annotation)[f], srt = 45
      )*1.5 - yspacing,
      label = if (f %in% as_ci) 
        paste("C.I.", colnames(x$tip.annotation)[1L])
      else
        colnames(x$tip.annotation)[f],
      pos = 1,
      srt = 45,
      xpd = NA
    )
    
    rect.args$ybottom <- min(tips[,2] - yspacing)
    rect.args$ytop    <- max(tips[,2] + yspacing)
    rect.args$col     <- "transparent"
    rect.args$border  <- "darkgray"
    rect.args$lwd     <- 1.5
    
    rect.args$xleft   <- (ndraws - 1)/nblocks + 1/nblocks*.05
    rect.args$xright  <- ndraws/nblocks - 1/nblocks*.05
    
    do.call(graphics::rect, rect.args)
    ndraws <- ndraws + 1L
    
    if (f %in% as_ci) {
      graphics::abline(v = (x1 - x0)/2 + x0, lty = 2)
      break
    }
    
  }
  
  # Drawing a legend
  width <- dev_size[1]*(1 - prop)
  graphics::par(op2)
  graphics::par(mai = c(0,0,dev_size[2] - op2$mai[1], dev_size[1]*prop + width/2))
  graphics::plot.window(c(0,1), c(0, 1))
  graphics::box(col="transparent")
  graphics::legend(
    "center",
    legend  = c("Duplication", "Other"),
    pt.bg   = node.type.col,
    pch     = 21,
    pt.cex  = 2,
    bty     = "n",
    horiz   = FALSE,
    title   = "Node type"
  )
  
  graphics::par(op2)
  graphics::par(mai = c(0, width/2, dev_size[2] - op2$mai[1], dev_size[1]*prop))
  graphics::plot.window(c(0,1), c(0,1))
  graphics::box(col="transparent")
  graphics::legend(
    "center",
    legend  = c("No function", "Function", "no information"),
    fill    = blue(c(0,1,.5)),
    bty     = "n",
    # density = c(NA, NA, 10),
    horiz   = FALSE,
    title   = "Annotations"
  )
  
  
  invisible(NULL)
}

