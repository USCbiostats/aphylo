#' Approximation of Geodesic distances using Matrix Powers
#' 
#' Given an adjacency matrix \eqn{A}, the geodesic can be approximated using
#' its powers, since each \eqn{(i,j)} element of \eqn{A^p} corresponds to the
#' number of \eqn{p} length steps between nodes \eqn{i} and \eqn{j}.
#' 
#' @template edges
#' @param nsteps Integer scalar. Number of maximum steps for the approximation.
#' @param warn Logical scalar. When \code{TRUE} shows a warning after no further
#' steps are needed.
#' @param undirected Logical scalar. When \code{TRUE} (default), the edgelist is treated
#' as undirected (see details).
#' @return A square matrix of size \eqn{n} with the shortest path between each
#' pair of nodes.
#' 
#' @details
#' When \code{undirected = TRUE}, the function extends \code{edges} such that
#' \code{edges = rbind(edges, edges[,2:1])}.
#' 
#' @author George G. Vega Yon
#' @references
#' This is a modified version of the function of the same name in the
#' R package \CRANpkg{netdiffuseR}.
#' @export
#' @name approx_geodesic
NULL


#' @export
#' @rdname approx_geodesic
approx_geodesic <- function(
  edges,
  nsteps     = 1e5,
  undirected = TRUE,
  warn       = FALSE
  ) 
  UseMethod("approx_geodesic")

#' @export
#' @rdname approx_geodesic
shortest_path <- approx_geodesic

#' @export
#' @rdname approx_geodesic
approx_geodesic.aphylo <- function(
  edges,
  nsteps     = 1e5,
  undirected = TRUE,
  warn       = FALSE
) {
  
  approx_geodesic.(edges$tree$edge - 1L, nsteps = nsteps, undirected = undirected,
                   warn = warn)
  
}

#' @export
#' @rdname approx_geodesic
approx_geodesic.phylo <- function(
  edges,
  nsteps     = 1e5,
  undirected = TRUE,
  warn       = FALSE
) {
  
  approx_geodesic.(edges$edge - 1L, nsteps = 1e5, undirected = undirected,
                   warn = warn)
  
}

#' @export
#' @rdname approx_geodesic
approx_geodesic.matrix <- function(
  edges,
  nsteps     = 1e5,
  undirected = TRUE,
  warn       = FALSE
) {
  
  if (min(edges) < 1)
    stop("`edges` must be integers starting from 1.", call. = FALSE)
  
  approx_geodesic.(edges - 1L, nsteps = 1e5, undirected = undirected,
                   warn = warn)
  
}
