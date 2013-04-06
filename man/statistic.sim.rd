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
statistic.sim(krige.obj, level, alternative = "less", ...)
}
\arguments{
  \item{krige.obj}{
An object from krige.uk.
}
  \item{level}{
The threshold/exceedance level under consideration.
}
  \item{alternative}{
Whether the region under consideration is for the values greater than level (alternative = "less"), for the values less than level (alternative = "greater"), or contour lines (alternative = "two.sided").  Defaults to "less".
}
  \item{...}{
Additional arguments when alternative = "two.sided".  See Details.
}
}

\details{
	When alternative = "two.sided", the "..." argument must include \code{user.cov} (a user-specified covariance function), \code{pgrid} (the grid of locations to be predicted, produced by \code{create.pgrid} or \code{create.pgrid2}), \code{X} (the matrix of covariates for the observed data), and any other arguments needed by \code{user.cov}.  Note that \code{user.cov} should take \code{cLcoords} as its first argument (a matrix containing the coordinates contour lines under consideration).  Additional arguments to \code{user.cov} are passed internally using the "..." argument.  The \code{user.cov} function should return a list with values V (the covariance matrix of the observed data), Vop (the cross-covariance matrix between the observed data and the responses with coordinates in cL), Vp (the covariance matrix of the responses with coordinates in cL), and Xp (the matrix of covariates for the coordinates contained in cL).  See below for an example.
}

\value{
A list containing the following objects:
  \item{statistic}{A vector with the observed values of the test statistic.}
  \item{statistic.sim}{A matrix containing the simulated values of the test statistic.}
  \item{alternative}{A character vector specifying the alternative hypothesis used to determine statistic.sim.}
  \item{level}{The threshold level used to calculate statistic and statistic.sim.}
  \item{level}{Additional arguments needed when alternative = "two.sided".  See Details.}
}

\author{
Joshua P. French
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
library(SpatialTools)	
	
set.seed(10)

# Load data
data(sdata)

# Create prediction grid
pgrid <- create.pgrid(0, 1, 0, 1, nx = 26, ny = 26)
pcoords <- pgrid$pgrid

# Create design matrices
X <- cbind(1, coords)
Xp <- cbind(1, pcoords)

# Generate covariance matrices V, Vp, Vop using appropriate parameters for observed
# data and responses to be predicted
spcov <- cov.sp(coords = coords, sp.type = "exponential", sp.par = c(1, 1.5), 
    error.var = 1/3, finescale.var = 0, pcoords = pcoords)

# Predict responses at pgrid locations
krige.obj <- krige.uk(y = as.vector(y), V = spcov$V, Vp = spcov$Vp, 
	Vop = spcov$Vop, X = X, Xp = Xp, nsim = 2000, 
	Ve.diag = rep(1/3, length(y)) , method = "chol")

# Simulate distribution of test statistic for different alternatives
statistic.sim.obj.less <- statistic.sim(krige.obj = krige.obj, 
	level = 5, alternative = "less")
statistic.sim.obj.greater <- statistic.sim(krige.obj = krige.obj, 
	level = 5, alternative = "greater")
	
# Construct null and rejection sets for two scenarios
n90 <- exceedance.ci(statistic.sim.obj.less, 
	conf.level = .90, type = "null")
r90 <- exceedance.ci(statistic.sim.obj.greater, 
	conf.level = .90, type = "rejection")

# Plot results
plot(pgrid, n90, col="blue", add = FALSE, xlab = "x", ylab = "y")
plot(pgrid, r90, col="orange", add = TRUE)
legend("bottomleft", 
	legend = c("contains true exceedance region with 90 percent confidence", 
	"is contained in true exceedance region with 90 percent confidence"), 
	col = c("blue", "orange"), lwd = 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ test statistic }

\references{
	French, J.P. and Sain, S.R. Spatio-Temporal Exceedance Locations and Confidence Regions. Annals of Applied Statistics.  Prepress.
}
