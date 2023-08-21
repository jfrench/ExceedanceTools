#' Construct confidence regions for exceedance (excursion) sets.
#' 
#' \code{confreg} constructs confidence regions for the exceedance (excursions) sets of geostatistical processes.  These will actually be credible regions if \code{obj} contains samples from the joint posterior predictive distribution in a Bayesian setting.
#' 
#' \code{obj} can be an object of class \code{matrix},  \code{krigeConditionalSample}, or \code{jointPredictiveSample}.  If \code{obj} is a \code{matrix}, then it should have \code{m} rows and \code{nsim} columns.  In that case, each row of \code{obj} corresponds to a sample from the conditional distribution of the response conditional on the observed data.  Each row represents a different location.  Generally, these locations are assumed to be on a grid spanning the spatial domain of interest.  A \code{krigeConditionalSample} object can be obtained using the \code{krige.sk}, \code{krige.ok}, 
#' or \code{krige.uk} functions in the \code{SpatialTools} package.  In these functions, the \code{nsim} argument must be greater than 0, and indicates the number of samples used to construct the confidence region.  A \code{jointPredictiveSample} object can be obtained using the \code{spLMPredictJoint} function in the \code{SpatialTools} package.  Since this is in the context of Bayesian statistics, the function actually produces credible region.  
#' 
#' If \code{statistic} is supplied for the direct construction procedure, then the locations are ordered by marginal probability and then the statistic.  \code{statistic} should be a vector of length \code{m}, where \code{m} is the number of prediction locations at which samples were drawn for in \code{obj}.
#' 
#' If type == \code{"o"}, then an outer credible region is constructed.  The outer credible region should entirely contain the true exceedanace region with the specified posterior probability.  If type == \code{"i"}, then an inner credible region is constructed.  The inner confidence region should be entirely contained within the true exceedanace region with specified posterior probability.
#' 
#' @param obj An object of the appropriate type (\code{matrix}, \code{krigeConditionalSample}, or \code{jointPredictiveSample}.  See Details.  
#' @param level The threshold level for the exceedance region.
#' @param statistic The statistic used in constructing the confidence region.  Should be a vector containing a value for each location
#' @param conf.level The confidence level of the confidence region.  Default is 0.95.
#' @param direction The direction of the exceedance region.  \code{">"} indicates the exceedance region is values above a threshold, while \code{"<"} indicates values below a threshold.
#' @param type \code{"o"} indicates on outer confidence region while \code{"i"} indicates in inner confidence region.
#' @param method \code{"test"} indicates a testing-based method, while \code{"direct"} indicates a direct method using joint probabilities.
#' @param greedy Only applicable for the direct construction method.  Default is \code{FALSE}.  If \code{TRUE}, then grid cells are added to the confidence region using a greedy algorithm based on joint probability. 
#' 
#' @return Returns an object of class \code{confreg} with the following components: 
#' \item{confidence}{The sites included in the confidence region.}
#' \item{complement}{The complement of the confidence region.}
#' 
#' @importFrom matrixStats colProds
#' @references 
#' Joshua P. French and Stephan R. Sain (2013).
#' Spatio-temporal exceedance locations and confidence
#' regions.  Annals of Applied Statistics.  7(3):1421-1449.
#'
#' French, J. P. (2014), Confidence regions for the level
#' curves of spatial data, Environmetrics, 25, pages
#' 498–512, DOI: 10.1002/env.2295
#'
#' French, J. P., and Hoeting, J. A. (2016) Credible regions
#' for exceedance sets of geostatistical data.
#' Environmetrics, 27: 4–14. doi: 10.1002/env.2371.
#' @author Joshua French
#' @export
#' @examples 
#' # Set parameters
#' n <- 100
#' mygrid = create.pgrid(0, 1, 0, 1, nx = 5, ny = 4)
#' n.samples <- 10
#' burnin.start <- 1
#' sigmasq <- 1
#' tausq <- 0.0
#' phi <- 1
#' cov.model <- "exponential"
#' n.report <- 5
#' 
#' # Generate coordinates
#' coords <- matrix(runif(2 * n), ncol = 2) 
#' pcoords <- mygrid$pgrid
#' # Construct design matrices
#' X <- as.matrix(cbind(1, coords))
#' Xp <- cbind(1, pcoords)
#' 
#' # Specify priors
#' starting <- list("phi" = phi, "sigma.sq"= sigmasq, "tau.sq" = tausq)
#' tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
#' priors.1 <- list("beta.Norm"=list(c(1, 2, 1), diag(100, 3)), "phi.Unif"=c(0.00001, 10), 
#'  "sigma.sq.IG"=c(1, 1))
#' 
#' # Generate data
#' library(SpatialTools)
#' B <- rnorm(3, c(1, 2, 1), sd = 10)
#' phi <- runif(1, 0, 10)
#' sigmasq <- 1/rgamma(1, 1, 1)
#' V <- simple.cov.sp(D = dist1(coords), cov.model, c(sigmasq, 1/phi), error.var = tausq, 
#'  smoothness = nu, finescale.var = 0)
#' y <- X %*% B + rmvnorm(1, rep(0, n), V) + rnorm(n, 0, sqrt(tausq))
#' 
#' # Create spLM object
#' library(spBayes)
#' m1 <- spBayes::spLM(y ~ X - 1, coords = coords, starting = starting, tuning = tuning, 
#'  priors = priors.1, cov.model = cov.model, n.samples = n.samples, verbose = FALSE,
#'  n.report = n.report)
#' 
#' # Sample from joint posterior predictive distribution
#' y1 <- spLMPredictJoint(m1, pred.coords = pcoords, pred.covars = Xp, 
#'  start = burnin.start, verbose = FALSE, method = "chol")
#'  u = quantile(y, .5)
#' myfun = function(x)
#' {
#'     (mean(x) - u)/sd(x)
#' }
#' 
#' myfun2 = function(x)
#' {
#'  mean(x > u)
#' }
#' 
#' stat1 = apply(y1, 1, myfun)
#' stat2 = apply(y1, 1, myfun2)
#' 
#' myconf = confreg(y1, level = u, statistic = NULL, direction = ">", type = "o", method = "direct")
#' myconf2 = confreg(y1, level = u, statistic = stat1, direction = ">", type = "o")
#' myconf3 = confreg(y1, level = u, statistic = stat2, direction = ">", type = "o")

confreg = function(obj, level, statistic = NULL, 
                   conf.level = 0.95, direction = ">", 
                   type = "o", method = "test", greedy = FALSE)
{
  if(!any(is.element(class(obj), c("matrix", "krigeConditionalSample", "jointPredictiveSample"))))
  {
    stop("obj is not of appropriate class")
  }
  if(is.element("krigeConditionalSample", class(obj))) obj = obj$sim
  m = nrow(obj)
  nsim = ncol(obj)
  if(length(level) > 1) stop("length(level)!=1")
  if(!is.numeric(level)) stop ("level must be a finite number")
  if(is.infinite(level)) stop("level must be a finite number")
  if(conf.level <= 0 | conf.level >=1) stop("conf.level must be between 0 and 1")
  if(!is.null(statistic))
  {
    statistic = c(statistic)
    if(length(statistic) != m) stop("length(statistic) != number of prediction locations")
  }
  if(!is.element(direction, c(">", "<")))
  {
    stop("Invalid direction")
  }
  if(!is.element(type, c("o", "i"))) stop("Invalid type")
  if(!is.element(method, c("test", "direct"))) stop("Invalid method")

  # inner confidence region is just the complement of the outer confidence
  # region for the other direction
  if(type == "i")
  {
    if(direction == ">"){ direction = "<" }
    else if(direction == "<"){ direction = ">" }
  }
  if(direction == "<"){ obj = -obj; level = -level }
  
  if(method == "test")
  {
    # determine quantile type (i.e., whether approximation is needed).
    # qtype <- 1 is exact
    qtype <- 2
    if((nsim * conf.level) == floor(nsim * conf.level)){ qtype <- 1 }

    statistic.sim = rep(Inf, nsim)
    
    for(j in 1:nsim)
    {
      which.exceedance <- which(obj[,j] >= level)
      if(length(which.exceedance) > 0)
      {
        statistic.sim[j] <- min(statistic[which.exceedance])
      }
    }
    cv <- stats::quantile(statistic.sim, prob = 1 - conf.level, type = qtype)
    confidence <- which(statistic >= cv)
    complement <- which(statistic < cv)
  }
  else if(method == "direct")
  {
    m = nrow(obj)
    if(direction == ">")
    {
      ind = obj < level
    }
    else if(direction == "<")
    {
      ind = obj > level
      if(!is.null(statistic)) statistic = -statistic
    }
    if(is.null(statistic)){
      ord <- order(rowMeans(ind), decreasing = TRUE)  
    }else
    {
      ord <- order(rowMeans(ind), statistic, decreasing = TRUE)  
    }
    
    if(!greedy)
    {
      A = apply(ind[ord,], 2, cumsum)
      ordind = ((apply(ind[ord,], 2, cumsum) - 
                   matrix(0:(m-1), nrow = m, ncol = nsim))
                == 1) + 0
      prob = rowMeans(ordind)
      complement = which(prob > conf.level)
      confidence = setdiff(1:m, complement)
      confidence = ord[confidence]
      complement = ord[complement]
    }else
    {
      complement = greedy.search(ind, conf.level)
      confidence = setdiff(1:m, complement)
    }

  }
  
  if(type == "i")
  {
    temp = confidence
    confidence = complement
    complement = temp
  }
  out = list(confidence = confidence, complement = complement)
  class(out) = "confreg"
  return(out)
}

greedy.search = function(ind, conf.level)
{
  rm_ind = rowMeans(ind)
  nsim = ncol(ind)
  potential = which(rm_ind >= conf.level)
  if(length(potential) == 0) return(numeric(0))
  always = which(rm_ind == 1)
  if(length(always) > 0){
    include = always
    potential = setdiff(potential, always)
  }else{
    include = which.max(rm_ind)
    potential = setdiff(potential, include)
  }
  
  continue = TRUE
  while(continue)
  {
    curProd = matrixStats::colProds(ind[include, , drop = FALSE]) # summarize current joint relationships
    newp = tcrossprod(curProd, ind[potential, , drop = FALSE])/nsim # calculate new joint probabilities
    if(max(newp) >= conf.level){
      newcell = potential[which.max(newp)]
      include = c(include, newcell)
      potential = setdiff(potential, newcell)
    }else
    {
      continue = FALSE
    }
  }
  include
}
