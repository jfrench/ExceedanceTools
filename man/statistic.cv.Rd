\name{statistic.cv}
\alias{statistic.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return critical value of distribution.
}
\description{
This function returns the critical value of the distribution of the test statistics based on the specified confidence level.  However, it is not recommended for general usage.  See note below.
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

\references{
French, J.P. and Sain, S.R. Spatio-Temporal Exceedance Locations and Confidence Regions. Annals of Applied Statistics.  Prepress.
}

\author{
Joshua P. French
}
\note{
This function can be used directly, but is used in the \code{exceedance.ci} function.  Using the exceedance.ci function to create confidence sets is the recommended approach.
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

# Calculate quantiles of distribution of statistic
q90.less <- statistic.cv(statistic.sim.obj.less, conf.level = .90)
q90.greater <- statistic.cv(statistic.sim.obj.greater, conf.level = .90)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cv}
\keyword{critical value}
\keyword{test statistic}

