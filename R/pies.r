pieslice <- function(a0, a1, r, d, x0, y0, s) {

  # Intermideate points
  mid <- seq(a0, a1, by = 2*pi/s*sign(a1 - a0))
  
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
  
  colnames(ans) <- c("x", "y")
  ans
  

}

circle <- function(x0, y0, r) {
  ans <-  pieslice(0, 2*pi, r=r, d=0, x0, y0, s=100)
  ans[-nrow(ans), ]
}


#' A flexible piechart.
#' 
#' @param x Numeric vector. Values that specify the area of the slices.
#' @param add Logical scalar. When \code{TRUE} it is added to the current device.
#' @param r Numeric scalar. Radious of the pie.
#' @param doughnut Numeric scalar. Radious of the inner circle (doughnut).
#' @param origin Numeric vector of length 2. Coordinates of the origin.
#' @param smooth Numeric scalar. Smoothness of the slices curve.
#' @param labs Character vector of length \code{length(x)}. Passed to
#' \code{\link[graphics:text]{text}}.
#' @param text.off Numeric scalar. Proportion that the text should be off (on)
#' the circle with respect to \code{r}.
#' @param text.args List. Further arguments passed to \code{\link[graphics:text]{text}}.
#' @param ... Further arguments passed to \code{\link[graphics:polygon]{polygon}}.
#' 
#' @return 
#' A list with two elements:
#' \item{slices}{A list of length \code{length(x)} with the coordinates of each
#'   slice.}
#' \item{textcoords}{A numeric matrix of size \code{length(x)*2} with 
#'   coordinates where the labels can be put at.}
#' 
#' @export
#' @examples
#'  
#' vals <- c(1,2,3,10)
#' piechart(vals, col=grDevices::blues9[5:8], border=NA, doughnut = .5, r=.75, labs=vals, text.off = .05)
#' piechart(vals, col=grDevices::blues9[3:6], border=NA, doughnut = .3, r=.5, add=TRUE)
#' piechart(vals, col=grDevices::blues9[1:4], border=NA, doughnut = .1, r=.3, add=TRUE)
piechart <- function(
  x,
  add = FALSE,
  r = 1,
  doughnut=0, 
  origin = c(0,0),
  smooth = 100,
  labs = NULL,
  text.off = .1,
  text.args = list(),
  ...) {
  
    # Assigning alpha
  alpha1 <- cumsum(x/sum(x)*2.0*pi)
  alpha0 <- c(0, alpha1[-length(x)])
  
  ans <- mapply(
    pieslice,
    a0 = alpha0,
    a1 = alpha1,
    r=r, d=doughnut, x0=origin[1], y0=origin[2],
    s = smooth, SIMPLIFY = FALSE
    )
  
  # Creating the device
  if (!add) {
    graphics::plot.new()
    graphics::plot.window(xlim=c(-r,r)*1.05, ylim = c(-r,r)*1.05)
  }
  
  # Adding the polygons
  mapply(graphics::polygon,
      x = lapply(ans, "[", j=1, i=),
      y = lapply(ans, "[", j=2, i=),
      ..., SIMPLIFY = FALSE
      )
  
  # Midpoints
  textcoords <- (alpha0 + alpha1)/2
  textcoords <- cbind(
    origin[1] + cos(textcoords)*r*(1 + text.off),
    origin[2] + sin(textcoords)*r*(1 + text.off)
    )
  
  # If labels are passed
  if (length(labs))
    do.call(
      graphics::text,
      c(
        list(x = textcoords[,1], y = textcoords[,2], labels=labs),
        text.args
        )
    )
  
  # Returning
  invisible(
      list(
        slices     = ans,
        textcoords = textcoords
        )
    )
  
}
  
  


