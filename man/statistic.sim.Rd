\name{statistic.sim}
\alias{statistic.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates distribution of test statistics for true exceedance region.
}
\description{
Simulates distribution of test statistics for true exceedance region.
}
\usage{
statistic.sim(krige.obj, level, alternative = "less")
}
\arguments{
  \item{krige.obj}{
An object from krige.uk.
}
  \item{level}{
The threshold/exceedance level under consideration.
}
  \item{alternative}{
Whether the region under consideration is for the values greater than level (alternative = "less") or for the values less than level (alternative = "greater").  Defaults to "less".
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
A list containing the following objects:
  \item{statistic}{A vector with the observed values of the test statistic.}
  \item{statistic.sim}{A matrix containing the simulated values of the test statistic.}
  \item{alternative}{A character vector specifying the alternative hypothesis used to determine statistic.sim.}
  \item{level}{The threshold level used to calculate statistic and statistic.sim.}
}

\author{
Joshua P. French
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
rnorm(10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
