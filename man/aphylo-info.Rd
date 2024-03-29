% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aphylo-ape-info.R
\name{aphylo-info}
\alias{aphylo-info}
\alias{Nann}
\alias{Nannotated}
\alias{Ntrees}
\title{Information about \code{aphylo} and \code{multiAphylo} objects}
\usage{
Nann(phy)

Nannotated(phy)

Ntrees(phy)
}
\arguments{
\item{phy}{Either an object of class \link{aphylo}, \link{multiAphylo}, or
\link{aphylo_estimates}.}
}
\value{
If \code{phy} is of class \code{aphylo}, then a single scalar.
otherwise, if \code{phy} is of class \code{multiAphylo}
}
\description{
Information about annotations, in particular, number of annotations (\code{Nann}),
number of annotated leaves (\code{Nannotated}), number of unnanotated leaves
(\code{Nunannotated}), and number of trees (\code{Ntrees}).
}
\examples{
# Generating data for the example
set.seed(223)
dat <- rmultiAphylo(10, n = 5, P = 2)
Nann(dat)
Nannotated(dat)
Ntrees(dat)
}
\seealso{
Other information: 
\code{\link{ape-methods}}
}
\concept{information}
