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

\references{
French, J.P. and Sain, S.R. Spatio-Temporal Exceedance Locations and Confidence Regions. Annals of Applied Statistics.  Prepress.
}

\author{
Joshua P. French
}

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
	Vop = spcov$Vop, X = X, Xp = Xp, nsim = 100, 
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
\keyword{exceedance}
\keyword{critical value}

