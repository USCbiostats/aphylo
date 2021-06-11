#' Rotation of polygon
#' @param mat Two-column numeric matrix. Coordinates of the polygon.
#' @param origin Numeric vector of length two. Origin.
#' @param alpha Numeric scalar. Rotation degree in radians.
#' @noRd
rotate <- function(mat, origin, alpha) {
  R <- matrix(
    c(cos(alpha), -sin(alpha), sin(alpha), cos(alpha)),
    nrow = 2, byrow = TRUE)

  origin <- matrix(c(origin[1], origin[2]), ncol=2, nrow = nrow(mat), byrow = TRUE)
  t(R %*% t(mat - origin)) + origin
}
