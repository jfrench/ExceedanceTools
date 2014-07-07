#' Construct confidence regions for exceedance (excurions) sets.
#' 
#' \code{confreg} constructs confidence regions for the exceedance (excursions) sets of geostatistical processes.
#' 
#' If type == \code{"o"}, then an outer confidence region is constructed.  The outer confidence region should entirely contain the true exceedanace region with high confidence.  If type == \code{"i"}, then an inner confidence region is constructed.  The inner confidence region should be entirely contained within the true exceedanace region with high confidence.
#' 
#' @param obj An object of class \code{jointPredictiveSample}.
#' @param level The threshold level for the exceedance region.
#' @param conf.level The confidence level of the confidence region.  Default is 0.95.
#' @param direction The direction of the exceedance region.  \code{">"} indicates the exceedance region is values above a threshold, while \code{"<"} indicates values below a threshold.
#' @param type \code{"o"} indicates on outer confidence region while \code{"i"} indicates in inner confidence region.
#' @param method If \code{method == 1}, then the region is constructed using a (test) statistic based procedure.  If \code{method == 2}, then the region is constructed directly.
#' @param statistic The type of statistic to use if \code{method == 1}.  If \code{statistic == 1}, then the statistic at each site is (prediction - level)/sqrt(mse).  If \code{statistic == 2}, then the statistic at each site is (prediction - level)/E[(Y - level)^2].  If \code{statistic == 3}, then the statistic is simply the predictive probability that the response is above (or below depending on context) the threshold \code{level}. 
#' @param ... Currently unimplemented.
#' 
#' @return Returns an object of class \code{confreg} with the following components: 
#' \item{confidence}{The sites included in the confidence region.}
#' \item{complement}{The complement of the confidence region.}
#' 
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
#' myconf = confreg(y1, level = quantile(y, .5), direction = ">", type = "o", method = 1, 
#'  statistic = 1)

confreg = function(obj, ...){
  UseMethod("confreg")
}

#' @rdname confreg
#' @export
confreg.jointPredictiveSample  = function(obj, level, conf.level = 0.95, 
                                          direction = ">", type = "o", 
                                          method = 1, statistic = 1,...)
{
  if(length(level) > 1) stop("length(level)!=1")
  if(!is.numeric(level)) stop ("level must be a finite number")
  if(is.infinite(level)) stop("level must be a finite number")
  if(conf.level <= 0 | conf.level >=1) stop("conf.level must be between 0 and 1")
  if(!is.element(direction, c(">", "<")))
  {
    stop("Invalid direction")
  }
  if(!is.element(type, c("o", "i"))) stop("Invalid type")
  if(!is.element(method, 1:2))
  {
    stop("method must equal 1 or 2")
  }
  if(!is.element(statistic, 1:3))
  {
    stop("method must equal 1, 2, or 3")
  }
  
  nsim = ncol(obj)
  
  # inner confidence region is just the complement of the outer confidence
  # region for the other direction
  if(type == "i")
  {
    if(direction == ">"){ direction = "<" }
    else if(direction == "<"){ direction = ">" }
    else{ stop("type == 'i' not valid for direction == '='") }
  }
  if(direction == "<"){ obj = -obj; level = -level }
  
  if(method == 1)
  {
    # determine quantile type (i.e., whether approximation is needed).
    # qtype <- 1 is exact
    qtype <- 2
    if((nsim * conf.level) == floor(nsim * conf.level)){ qtype <- 1 }
    
    # if a sample doesn't have an exceedance, the extreme statistic
    # for the sample should be infinite for statistics 1 and 2.
    # it should be 0 for statistic 3.
    if(statistic == 1 | statistic == 2) statistic.sim = rep(-Inf, nsim)
    else{ statistic.sim = rep(0, nsim) }
    
    pred = apply(obj, 1, mean)
    if(statistic == 1){ stat = (rowMeans(obj) - level)/apply(obj, 1, sd) }
    else if(statistic == 3){ stat = apply(obj >= level, 1, mean) }
    else{ stat = (rowMeans(obj) - level)/sqrt(rowMeans((obj - level)^2)) }
    
    for(j in 1:nsim)
    {
      which.exceedance <- which(obj[,j] >= level)
      if(length(which.exceedance) > 0)
      {
        statistic.sim[j] <- min(stat[which.exceedance])
      }
    }
    cv <- quantile(statistic.sim, prob = 1 - conf.level, type = qtype)
    confidence <- which(stat >= cv)
    complement <- which(stat < cv)
  }
  else if(method == 2)
  {
    m = nrow(obj)
    ind = obj < level
    ord <- order(rowMeans(ind), decreasing = TRUE)
    A = apply(ind[ord,], 2, cumsum)
    ordind = ((apply(ind[ord,], 2, cumsum) - 
                 matrix(0:(m-1), nrow = m, ncol = nsim))
              == 1) + 0
    prob = rowMeans(ordind)
    complement = which(prob > conf.level)
    confidence = setdiff(1:m, complement)
    confidence = ord[confidence]
    complement = ord[complement]
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
