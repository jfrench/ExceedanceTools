\name{exceedance.ci}
\alias{exceedance.ci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return confidence set.
}
\description{
Returns a confidence set for an exceedance region.
}
\usage{
exceedance.ci(statistic.sim.obj, conf.level = 0.95, type = "null")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{statistic.sim.obj}{
An object from the \code{statistic.sim} object.
}
  \item{conf.level}{
The desired confidence level of the confidence region.
}
  \item{type}{
A character vector of 'null' or 'rejection' indicating whether the function should return the null region or rejection region, respectively, related to \code{statistic.sim.obj}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a numeric vector with the set of pixels comprising the null or rejection region related to statistic.sim.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Joshua P. French
}

\examples{
rnorm(10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{exceedance}
\keyword{critical value}

