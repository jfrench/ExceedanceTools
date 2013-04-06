\name{create.pgrid2}
\alias{create.pgrid2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create grid of locations using vectors of x- and y-coordiantes}
\description{
Create grid of locations using vector of x-coordinates and y-coordinates.
}
\usage{
create.pgrid2(xgrid, ygrid, midpoints = FALSE, poly.coords = NULL)
}
\arguments{
  \item{xgrid}{
A vector of locations in the x direction.
}
  \item{ygrid}{
A vector of locations in the y direction.
}
  \item{midpoints}{
A boolean value (TRUE or FALSE) indicating whether the boundary values are for the midpoint of a pixel (midpoints = TRUE) or for the boundary of the spatial domain in general (midpoints = FALSE, in which case the midpoints are calculated internally).  Default is FALSE.}
  \item{poly.coords}{
An \eqn{n \times 2} matrix with the coordinates specifying the polygon vertices of the true spatial domain of interest within the rectangular boundaries provided by xmin, xmax, ymin, and ymax.  If this is provided, the grid returned will be within the convex hull of poly.coords.
}
}
\details{
The key argument in the function midpoints.  If this is TRUE, it is assumed that the boundaries of the spatial domain correspond to the midpoints of the cell/pixel in the grid.  Otherwise, it is assumed that the boundaries correspond to the actual borders of the region of interest.

If poly.coords is supplied, the grid returned is the grid of midpoints contained in the convex hull of poly.coords.}
\value{
  \item{pgrid}{An \eqn{np \times 2} matrix of locations (the midpoints of the pixelized grid).}
  \item{np}{The number of rows in pgrid.}
  \item{p.in.grid}{A vector of 0s and 1s indicating whether the midpoint of each pixel is in the convex hull of poly.coords.  If poly.coords is not provided, this is a vector of 1s.}
  \item{ubx}{The pixel boundaries in the x-direction.}
  \item{uby}{The pixel boundaries in the y-direction.}
  \item{upx}{The pixel midpoints in the x-direction.}
  \item{upy}{The pixel midpoints in the y-direction.}
}

\author{
Joshua French
}

\seealso{
pgrid
}
\examples{

pgrid <- create.pgrid2(seq(0, 1, len = 101), seq(0, 1, len = 101), midpoint = FALSE)
pgridb <- create.pgrid2(seq(.005, .995, len = 100), seq(.005, .995, len = 100), midpoint = TRUE)

# pgrids produced match
range(pgrid$pgrid - pgridb$pgrid)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{grid}
