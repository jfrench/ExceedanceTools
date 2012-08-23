\name{plot.pgrid}
\alias{plot.pgrid}
\title{
Plot grid of locations.
}
\description{
Plot grid of locations based on grid.  
}
\usage{

\method{plot}{pgrid}(x, set, col = "gray", add = FALSE, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The pgrid object (created by \code{create.pgrid})  from which you want to plot certain pixels.
}
  \item{set}{
Indicates the pixels/cells that should be plotted.
}
  \item{col}{
The color of the plotted pixels.
}
  \item{add}{
A boolean value indicating whether the pixels should be added to an existing plot (add = TRUE) or should the pixels be plotted on a new plot (add = FALSE).
}
  \item{\dots}{
Additional arguments that will be passed to the plot function (assuming add = FALSE).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
This function does not return anything; it only creates a new plot or modifies an existing plot.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Joshua French
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
rnorm(10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}

