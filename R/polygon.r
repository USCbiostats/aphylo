#' n-sided polygons
#' Calculate the coordinates for an nsided polygon
#' @param x,y Numeric scalar. Origin of the polygon.
#' @param n Integer scalar. Number of sides.
#' @param r Numeric scalar. Radious of the polygon.
#' @param d Numeric scalar. Starting degree in radians.
#' @return
#' A two column matrix with the coordinates to draw a n sided polygon.
#' @examples
#' graphics.off()
#' oldpar <- par(no.readonly = TRUE)
#'
#' par(xpd = NA, mfrow = c(3, 3), mai = rep(0, 4))
#' for (n in c(2, 3, 4, 5, 6, 8, 12, 20, 50)) {
#'
#'   plot.new()
#'   plot.window(c(-1.25,1.25), c(-1.25,1.25))
#'
#'   for (i in seq(1, .0005, length.out = 200)) {
#'     col <- adjustcolor("tomato", alpha.f = i)
#'     polygon(npolygon(x=(i-1)/4, y = (i-1)/4, r = i, d = i-1, n = n),
#'             col = NA, border=col)
#'   }
#'
#'   mtext(sprintf("n = %i", n), side = 1, line = -3)
#' }
#'
#' par(oldpar)
#' @noRd
npolygon <- function(
  x = 0,
  y = 0,
  n = 6L,
  r = 1.0,
  d = 2.0*pi/(n)/2
) {

  deg <- seq(d, 2.0*pi + d, length.out = n + 1)[-(n+1)]

  coords <- cbind(
    x = cos(deg)*r + x,
    y = sin(deg)*r + y
  )

}
