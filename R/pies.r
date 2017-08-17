pieslice <- function(a0, a1, r, d, x0, y0, edges, off) {

  # Intermideate points
  mid <- seq(a0, a1, by = 2*pi/edges*sign(a1 - a0))
  
  # In case that the points are not sufficient
  if (length(mid) < 3)
    mid <- seq(a0, a1, length.out = 3)
  if (utils::tail(mid, 1) != a1)
    mid <- seq(a0, a1, length.out = length(mid))
  
  # Computing midpoints
  mid <- cbind(cos(mid), sin(mid))
  m <- nrow(mid)

  if (d != 0) {
    pbase <- (mid*d)[m:1,]
  } else pbase <- c(0,0)
  
  ans <- rbind.data.frame(
    mid*r, 
    pbase,
    make.row.names=FALSE
  ) 
  
  # Translating to the origin
  ans[,1] <- x0 + ans[,1]
  ans[,2] <- y0 + ans[,2]
  
  # Setting off
  ans[,1] <- ans[,1] + cos((a1 + a0)/2)*off
  ans[,2] <- ans[,2] + sin((a1 + a0)/2)*off
  
  colnames(ans) <- c("x", "y")
  ans
  

}

circle <- function(x0, y0, r) {
  ans <-  pieslice(0, 2*pi, r=r, d=0, x0, y0, edges=100, off=0)
  ans[-nrow(ans), ]
}

#' A flexible piechart.
#' 
#' @param x Numeric vector. Values that specify the area of the slices.
#' @param add Logical scalar. When \code{TRUE} it is added to the current device.
#' @param radius Numeric scalar. Radious of the pie.
#' @param doughnut Numeric scalar. Radious of the inner circle (doughnut).
#' @param origin Numeric vector of length 2. Coordinates of the origin.
#' @param edges Numeric scalar. Smoothness of the slices curve.
#' @param slice.off Numeric vector. When \code{!=0}, specifies how much to
#' move the slice away from the origin. When scalar is recycled.
#' @param labels Character vector of length \code{length(x)}. Passed to
#' \code{\link[graphics:text]{text}}.
#' @param tick.len Numeric scalar. Size of the tick marks as \% of the radious.
#' @param text.args List. Further arguments passed to \code{\link[graphics:text]{text}}.
#' @param segments.args List. Further arguments passed to \code{\link[graphics:segments]{segments}}
#' when drawing the tickmarks.
#' @param init.angle Numeric scalar. Angle from where to start drawing in degrees.
#' @param last.angle Numeric scalar. Angle where to finish drawing in degrees.
#' @param ... Further arguments passed to \code{\link[graphics:polygon]{polygon}}.
#' 
#' @return 
#' A list with two elements:
#' \item{slices}{A list of length \code{length(x)} with the coordinates of each
#'   slice.}
#' \item{textcoords}{A numeric matrix of size \code{length(x)*2} with 
#'   coordinates where the labels can be put at.}
#' \item{alpha0}{A numeric vector of size \code{length(x)} with the starting
#' degree in radians of the slice.}
#' \item{alpha1}{A numeric vector of size \code{length(x)} with the ending
#' degree in radians of the slice.}
#' 
#' @export
#' @examples
#'  
#' vals <- c(1,2,3,10)
#' piechart(vals, col=grDevices::blues9[5:8], border=NA, doughnut = .5, radius=.75, labs=vals)
#' piechart(vals, col=grDevices::blues9[3:6], border=NA, doughnut = .3, radius=.5, add=TRUE)
#' piechart(vals, col=grDevices::blues9[1:4], border=NA, doughnut = .1, radius=.3, add=TRUE)
piechart <- function(
  x,
  labels = names(x),
  radius = 1,
  doughnut = 0, 
  origin = c(0,0),
  edges = 200,
  slice.off = 0,
  init.angle = 0,
  last.angle = 360,
  tick.len = .1,
  text.args = list(),
  segments.args = list(),
  skip.plot.slices=FALSE,
  add = FALSE,
  ...) {
  
  # Assigning alpha
  init.angle <- init.angle/360*2.0*pi # as radians
  last.angle <- last.angle/360*2.0*pi
  
  alpha1 <- cumsum(x/sum(x)*last.angle) + init.angle
  alpha0 <- c(0, alpha1[-length(x)])
  
  ans <- mapply(
    pieslice,
    a0 = alpha0,
    a1 = alpha1,
    r = radius, d=doughnut, x0=origin[1], y0=origin[2],
    edges = edges, 
    off = slice.off, 
    SIMPLIFY = FALSE
    )
  
  # Fetching size
  cex <- if (length(labels) && ("cex" %in% names(text.args)))
    text.args[["cex"]]
  else 1
  
  # Creating the device
  if (!add) {
    graphics::plot.new()
    
    # Adjusting so that we get nice circles
    adj <- graphics::par()$pin
    adj <- adj[1]/adj[2]
    
    # Adjusting the size... including the labels
    ran <- if (length(labels))
      max(
        c(
          graphics::strwidth(labels, units="user", cex=cex),
          graphics::strheight(labels, units="user", cex=cex)
          ),
        na.rm=TRUE)*2
    else 0
    
    ran <- (ran + radius*1.1 + tick.len/2 + max(slice.off))*c(-1,1)

    if (adj > 1)
      graphics::plot.window(xlim=ran*adj, ylim = ran)
    else
      graphics::plot.window(xlim=ran, ylim = ran/adj)
  }
  
  # Adding the slices
  if (!skip.plot.slices)
    mapply(graphics::polygon,
        x = lapply(ans, "[", j=1, i=),
        y = lapply(ans, "[", j=2, i=),
        ..., SIMPLIFY = FALSE
        )
    
  # Midpoints
  angles     <- (alpha0 + alpha1)/2
  textcoords <- cbind(
    origin[1] + cos(angles)*(radius*1.05 + tick.len/2 + slice.off) ,
    origin[2] + sin(angles)*(radius*1.05 + tick.len/2 + slice.off)
    )
  
  # If labels are passed
  if (length(labels)) {
    
    # Adjusting according to string lenght
    textcoords <- textcoords +
      cbind(
        cos(angles)*graphics::strwidth(labels, cex=cex)/2,
        sin(angles)*graphics::strwidth(labels, cex=cex)/2
      )
    
    # Drawing the text
    do.call(mapply,
            c(list(FUN = graphics::text,
                 x=textcoords[,1],
                 y=textcoords[,2],
                 labels=labels,SIMPLIFY = FALSE),
                 text.args
            ))
    
    
    # Here should go the tick marks...
    x0 <- origin[1] + cos(angles)*radius*(1 - tick.len/2 + slice.off)
    x1 <- origin[1] + cos(angles)*radius*(1 + tick.len/2 + slice.off)
    y0 <- origin[2] + sin(angles)*radius*(1 - tick.len/2 + slice.off)
    y1 <- origin[2] + sin(angles)*radius*(1 + tick.len/2 + slice.off)
    
    toplot <- which(!is.na(labels))
    do.call(graphics::segments, 
            c(
              list(
                x0 = x0[toplot],
                y0 = y0[toplot],
                x1 = x1[toplot],
                y1 = y1[toplot]
              ),
              segments.args
    ))
    
  }
    
  
  # Returning
  invisible(
      list(
        slices     = ans,
        textcoords = textcoords,
        alpha0     = alpha0,
        alpha1     = alpha1
        )
    )
  
}
  

colorkey <- function(
  x0,y0,x1,y1,
  label.from,
  label.to,
  tick.range = c(0, 1),
  tick.marks = c(0, .25, .5, .75, 1),
  cols = c("white", "steelblue"),
  nlevels = 100
  ) {
  
  # Adjusting to textsize
  x0 <- x0 + graphics::strwidth(label.from)
  x1 <- x1 - graphics::strwidth(label.to)
  
  # Writing labels
  graphics::text(
    x = c(
      x0 - graphics::strwidth(label.from)/2, 
      x1 + graphics::strwidth(label.to)/2 
    ),
    y = rep((y1+y0)/2, 2),
    labels = c(label.from, label.to)
  )
  
  # Readjusting for more space
  x0 <- x0 + (x1-x0)/40
  x1 <- x1 - (x1-x0)/40
  
  # Computing coordinates
  cols <- grDevices::colorRampPalette(cols)(nlevels)
  n  <- length(cols)
  xsize <- (x1 - x0)/n
  
  xcoords <- seq(x0 + xsize/2, x1 - xsize/2, length.out = n)
  ycoords <- rep((y1+y0)/2, n)
  
  # xcoords[1] <- xcoords[1] + xsize
  
  # Drawing rectangles
  graphics::symbols(
    x = xcoords,
    y = ycoords,
    inches=FALSE,
    bg = cols,
    fg = "transparent",
    rectangles = cbind(rep(xsize, n), y1-y0),
    add=TRUE
  )
  
  # Bouding box
  graphics::symbols(
    x = (x1 + x0)/2,
    y = (y1 + y0)/2,
    rectangles = cbind((x1 - x0), y1 - y0),
    add=TRUE, inches = FALSE
  )
  
  # Top tickmarks
  tick.pos <- (tick.marks - tick.range[1])/(tick.range[2] - tick.range[1])*(x1-x0) + x0
  graphics::segments(
    x0 = tick.pos,
    y0 = y1 - (y1 - y0)/5 ,
    x1 = tick.pos,
    y1 = y1 + (y1 - y0)/5
  )
  
  graphics::text(
    x = tick.pos,
    y = y1 + (y1 - y0)/5 + graphics::strheight(max(tick.marks, na.rm=TRUE))/1.5,
    labels = tick.marks
  )
  
}

