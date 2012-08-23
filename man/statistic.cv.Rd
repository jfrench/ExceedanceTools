\name{statistic.cv}
\alias{statistic.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return critical value of distribution.
}
\description{
This function returns the critical value of the distribution of the test statistics based on the specified confidence level.
}
\usage{
statistic.cv(statistic.sim.obj, conf.level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{statistic.sim.obj}{
An object returned from the \code{statistic.sim} function.
}
  \item{conf.level}{
The desired confidence level of the confidence interval we want to construct.
}
}
%\details{
%}
\value{
Returns the desired critical value.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Joshua P. French
}
\note{
This function can be used directly, but is used in the \code{exceedance.ci} function.  Using the exceedance.ci function is the recommended approach.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
rnorm(10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cv}
\keyword{critical value}
\keyword{test statistic}

