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
\keyword{plot}

