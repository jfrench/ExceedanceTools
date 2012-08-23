\name{create.pgrid}
\alias{create.pgrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create grid of locations
}
\description{
Create grid of locations from boundaries of domain and other information.
}
\usage{
create.pgrid(xmin, xmax, ymin, ymax, nx, ny, midpoints = FALSE, poly.coords = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xmin}{
The minimum value of the boundary of the x-coordinates of the spatial domain.
}
  \item{xmax}{
The maximum value of the boundary of the x-coordinates of the spatial domain.
}
  \item{ymin}{
The minimum value of the boundary of the y-coordinates of the spatial domain.
}
  \item{ymax}{
The maximum value of the boundary of the y-coordinates of the spatial domain.
}
  \item{nx}{
The number of gridpoints/cells/pixels in the x direction.
}
  \item{ny}{
The number of gridpoints/cells/pixels in the y direction.
}
  \item{midpoints}{
A boolean value (TRUE or FALSE) indicating whether the boundary values are for the midpoint of a pixel (midpoints = TRUE) or for the boundary of the spatial domain in general (midpoints = FALSE, in which case the midpoints are calculated internally).  Default is FALSE.
}
  \item{poly.coords}{
An \eqn{n \times 2} matrix with the coordinates specifying the polygon vertices of the true spatial domain of interest within the rectangular boundaries provided by xmin, xmax, ymin, and ymax.  If this is provided, the grid returned will be within the convex hull of poly.coords.
}
}
\details{
The key argument in the function midpoints.  If this is TRUE, it is assumed that the boundaries of the spatial domain correspond to the midpoints of the cell/pixel in the grid.  Otherwise, it is assumed that the boundaries correspond to the actual borders of the region of interest.

If poly.coords is supplied, the grid returned is the grid of midpoints contained in the convex hull of poly.coords.
}
\value{
  \item{pgrid}{An \eqn{np \times 2} matrix of locations (the midpoints of the pixelized grid).}
  \item{np}{The number of rows in pgrid.}
  \item{p.in.grid}{A vector of 0s and 1s indicating whether the midpoint of each pixel is in the convex hull of poly.coords.  If poly.coords is provided, this is a vector of 1s.}
  \item{ubx}{The pixel boundaries in the x-direction.}
  \item{uby}{The pixel boundaries in the y-direction.}
  \item{upx}{The pixel midpoints in the x-direction.}
  \item{upy}{The pixel midpoints in the y-direction.}
}

\author{
Joshua French
}

\seealso{
pgrid2
}
\examples{
x.min <- y.min <- 0
x.max <- y.max <- 1
ngrid <- 50


pgrid <- create.pgrid(x.min, x.max, y.min, y.max, nx = ngrid, ny = ngrid, midpoint = FALSE)

pgridb <- create.pgrid(.01, .99, .01, .99, nx = ngrid, ny = ngrid, midpoint = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{grid}
