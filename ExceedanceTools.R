library(SpatialTools)

cov.condnorm <- function (y, V, Vp, Vop, method = "eigen", return.decomp = FALSE) 
{
    ViVop <- solve(V, Vop)
    Vc <- Vp - crossprod(Vop, ViVop)
    if(return.decomp)
    {
    	decomp.Vc <- decomp.cov(Vc, method = method)
    	return(decomp.Vc)
    }else
    {
	    return(Vc)	    
    }
}

ekrige.uk <- function(y, V, Vp, Vop, X, Xp)
{
	###compute matrix products for future use
	#compute Vi*X
	ViX <- solve(V, X)
	#compute X'*Vi*X
	XtViX <- crossprod(ViX, X)
	
	#compute gls estimates of regression coefficients
	coeff <- solve(XtViX, crossprod(ViX, y))

	ViVop <- solve(V, Vop)

	#compute kriging weights
	w <- solve(V, Vop - X%*%solve(XtViX, crossprod(X, ViVop) - t(Xp)))

	#blup for Yp
	pred <- crossprod(w, y)
	
	#variance of (Yp - pred)
	mspe <- colSums((V %*% w) * w) - 2 * colSums(w * Vop) + diag(Vp)

	out <- list(pred = pred, mspe = mspe, coeff = coeff, 
		vcov.coeff = solve(XtViX), ViVop = ViVop)

	return(out)
}

statistic.sim <- function(nsim = 1, y, X, Xp, ekrige.obj, decomp.Vc, level, alternative = "less", 
	return.yp.sim = FALSE)
{
	y <- as.vector(y)

	n <- length(y)
	nk <- length(ekrige.obj$coeff)
	np <- length(ekrige.obj$pred)
	statistic.sim <- numeric(nsim)
	
	coeff <- ekrige.obj$coeff
	decomp.vcov.coeff <- decomp.cov(ekrige.obj$vcov.coeff, method = "chol")
	
	coeff.sim <- matrix(coeff, nrow = nk, ncol = nsim) + 
		sqrt(2) * decomp.vcov.coeff %*% matrix(rnorm(nk * nsim), nrow = nk, ncol = nsim)
	ymat <- matrix(y, nrow = length(y), ncol = nsim)

	#conditional mean of unobserved Y given observed Y = y
	mcmat <- Xp %*% coeff.sim + crossprod(ekrige.obj$ViVop, ymat - X %*% coeff.sim)
	
	yp.sim <- mcmat + decomp.Vc %*% matrix(rnorm(nsim * np), nrow = np, ncol = nsim)

	statistic <- (ekrige.obj$pred - level)/sqrt(ekrige.obj$mspe)

	if(alternative == "less")
	{

		for(i in 1:nsim)
		{
			which.exceedance.sim <- which(yp.sim[, i] >= level)
			if(length(which.exceedance.sim) > 0)
			{
				statistic.sim[i] <- min(statistic[which.exceedance.sim])
			}
		}
	}else if(alternative == "greater")
	{
		for(i in 1:nsim)
		{
			which.exceedance.sim <- which(yp.sim[, i] <= level)
			if(length(which.exceedance.sim) > 0)
			{
				statistic.sim[i] <- max(statistic[which.exceedance.sim])
			}
		}
	}

	if(return.yp.sim)
	{
		return(list(statistic = statistic, statistic.sim = statistic.sim, alternative = alternative, 
			level = level, yp.sim = yp.sim))
	}else
	{
		return(list(statistic = statistic, statistic.sim = statistic.sim, alternative = alternative, 
			level = level))
	}
}

val2 <- function(i, u = NULL)
{
	rall <- ma + decomp.Va %*% rnorm(nrow(ma))
	
	y <- rall[opos,] + rnorm(length(opos), sd = sqrt(error.var))
	yp <- rall[ppos,]

	if(is.null(u)) { u <- quantile(yp, prob = lev/100) }

	true.exceedance <- which(yp > u)

	pred <- crossprod(w, y)
	statistic <- (pred - u)/sqrt(mspe)	
	
	#compute gls estimates of regression coefficients
	coeff <- solve(XtViX, crossprod(ViX, y))
	
	ekrige.obj <- list(pred = pred, mspe = mspe, coeff = coeff, vcov.coeff = vcov.coeff, 
		ViVop = ViVop)

	statistic.sim.obj <- statistic.sim(nsim = nsim, y = y, X = X, Xp = Xp, 
		ekrige.obj = ekrige.obj, decomp.Vc = decomp.Vc, level = u, 
		alternative = "less", return.yp.sim = FALSE)
	
	#Critical values of statistic
	q90 <- statistic.cv(statistic.sim.obj, conf.level = .90)
	q95 <- statistic.cv(statistic.sim.obj, conf.level = .95)

	#Construct null and rejection sets 
	n90 <- exceedance.ci(statistic.sim.obj, conf.level = .90, type = "null")
	n95 <- exceedance.ci(statistic.sim.obj, conf.level = .95, type = "null")	
	r90 <- exceedance.ci(statistic.sim.obj, conf.level = .90, type = "rejection")
	r95 <- exceedance.ci(statistic.sim.obj, conf.level = .95, type = "rejection")	

	#Determine whether rejection region intersects true exceedance region
	conf.success90 <- exceedance.success(r90, true.exceedance)
	conf.success95 <- exceedance.success(r95, true.exceedance)

	#Is the simulation a success
	success90 <- conf.success90$success
	success95 <- conf.success95$success
	
	#Proportion of true exceedance region contained in true confidence region
	pcontain90 <- conf.success90$per
	pcontain95 <- conf.success95$per

	if(i %% nreport==0){ print(i); flush.console()}

	return(list(statistic = statistic, true.exceedance = true.exceedance, 
		q90 = q90, q95 = q95, n90 = n90, n95 = n95,
		r90 = r90, r95 = r95, success90 = success90, success95 = success95,
		pcontain90 = pcontain90, pcontain95 = pcontain95, u = u, statistic.sim.obj))
}

#Pixelize the domain into n x n equal sized pixels
pixel.domain=function(xmin, xmax, ymin, ymax ,n, plot=FALSE)
{
	xstep=(xmax-xmin)/n	#Calculate the pixel width
	ystep=(ymax-ymin)/n	#Calculate the pixel height
	half.xstep=xstep/2	#Calculate half the pixel width
	half.ystep=ystep/2	#Calcualte half the pixel heiht

  mid.x=xmin+half.xstep+0:(n-1)*xstep #Calculates the x coordinate of the
                                      #centerpoint of each pixel
  mid.y=ymin+half.ystep+0:(n-1)*ystep #Calculates the y coordinate of the
                                      #centerpoint of each pixel
  mid.grid=expand.grid(mid.x,mid.y, KEEP.OUT.ATTRS = TRUE) #Make a grid of the centerpoint
                                    #of all of the pixels
  x.edges=xmin+0:n*xstep  #The x-coords of the pixels
  y.edges=ymin+0:n*ystep  #The y-coords of the pixels
  #grid=expand.grid(x.edges,y.edges)  #Makes a grid of the endpoints
                                     #of each pixel
  #Make a matrix of the x and y coordinates of the vertices of each pixel
  vert.coords=cbind(mid.grid[,1]-half.xstep,mid.grid[,1]+half.xstep,
                    mid.grid[,2]-half.ystep,mid.grid[,2]+half.ystep)

  out=NULL
  out$mid.grid=mid.grid
  out$x.edges=x.edges
  out$y.edges=y.edges
  out$vert.coords=vert.coords
  out
}

create.pgrid <- function(xmin, xmax, ymin, ymax, nx, ny, midpoints = TRUE,
	coords = NULL)
{
	if(midpoints)
	{
		xstep <- (xmax-xmin)/nx	#Calculate the pixel width
		ystep <- (ymax-ymin)/ny	#Calculate the pixel height
	}
	#half.xstep=xstep/2	#Calculate half the pixel width
	#half.ystep=ystep/2	#Calcualte half the pixel heiht

	#Determine x and y midpoints of all of the pixels	
	upx <- xmin + xstep/2 + 0:(nx-1) * xstep 
	upy <- ymin + ystep/2 + 0:(ny-1) * ystep 
	
	#Create boundaries for pixels
	ubx <- xmin + 0:nx * xstep
	uby <- ymin + 0:ny * ystep
	
	#If coords are supplied, create pgrid based on whether points
	#are contained in convex hull of coords.
	if(!is.null(coords))
	{
		#Determine points bounding convex hull of observed coordinates
		convex.hull <- chull(coords)
		#Create convex border
		convex.border <- coords[convex.hull, ]
		#Determine points of rectangular grid (based on xgrid and ygrid)
		#within convex border of observed coordinates
		convex.polygrid <- polygrid(upx, upy, convex.border, 
			vec.inout = TRUE)
		#Extract prediction coordinates within border
		pcoords <- convex.polygrid$xypoly
		#Determine number of prediction locations
		np <- nrow(pcoords)
		#Determine which points are a prediction coordinates within 
		#rectangular grid
		p.in.grid <- convex.polygrid$vec.inout
	}else
	{
		pcoords <- as.matrix(expand.grid(upx, upy))
		np <- length(upx) * length(upy)
		p.in.grid <- rep(TRUE, np)
	}
	
	return(list(pcoords = as.matrix(pcoords), np = np, p.in.grid = p.in.grid, 
		ubx = ubx, uby = uby, upx = upx, upy = upy))	
}


create.pgrid2 <- function(xgrid, ygrid, midpoints = FALSE, coords = NULL)
{
	nx <- length(xgrid)
	ny <- length(ygrid)
	
	if(!midpoints)
	{
		#Rename xgrid and ygrid (for later consistency).  
		ubx <- xgrid; uby <- ygrid

		#Determine x and y midpoints of all of the pixels
		upx <- (xgrid[-nx] + xgrid[-1])/2
		upy <- (ygrid[-ny] + ygrid[-1])/2
	}else
	{
		#Rename xgrid and ygrid (for later consistency)
		upx <- xgrid; upy <- ygrid

		#Create boundaries for pixels
		ubx <- (xgrid[-nx] + xgrid[-1])/2
		ubx <- c(2 * upx[1] - ubx[1], ubx, 2 * upx[nx] - ubx[nx - 1])

		uby <- (ygrid[-ny] + ygrid[-1])/2
		uby <- c(2 * upy[1] - uby[1], uby, 2 * upy[ny] - uby[ny - 1])
	}
	
	#If coords are supplied, create pgrid based on whether points
	#are contained in convex hull of coords.
	if(!is.null(coords))
	{
		#Determine points bounding convex hull of observed coordinates
		convex.hull <- chull(coords)
		#Create convex border
		convex.border <- coords[convex.hull, ]
		#Determine points of rectangular grid (based on xgrid and ygrid)
		#within convex border of observed coordinates
		convex.polygrid <- polygrid(upx, upy, convex.border, 
			vec.inout = TRUE)
		#Extract prediction coordinates within border
		pcoords <- convex.polygrid$xypoly
		#Determine number of prediction locations
		np <- nrow(pcoords)
		#Determine which points are a prediction coordinates within 
		#rectangular grid
		p.in.grid <- convex.polygrid$vec.inout
	}else
	{
		pcoords <- as.matrix(expand.grid(upx, upy))
		np <- length(upx) * length(upy)
		p.in.grid <- rep(TRUE, np)
	}
	
	return(list(pcoords = as.matrix(pcoords), np = np, p.in.grid = p.in.grid, 
		ubx = ubx, uby = uby, upx = upx, upy = upy))	
}

plot.pgrid <- function(set, pgrid, col = "gray", add = FALSE, ...)
{
	if(is.null(pgrid$ubx) || is.null(pgrid$uby))
	{
		stop("pgrid must contain boundaries.  Use return.boundaries = TRUE when running create.pgrid() function")
	}
	nx <- length(pgrid$ubx)
	ny <- length(pgrid$uby)
	
	#grid point bound for nw, ne, se, and sw of pcoords
	nw <- as.matrix(expand.grid(pgrid$ubx[-1], pgrid$uby[-ny]))[pgrid$p.in.grid,]
	ne <- as.matrix(expand.grid(pgrid$ubx[-nx], pgrid$uby[-ny]))[pgrid$p.in.grid,]
	se <- as.matrix(expand.grid(pgrid$ubx[-nx], pgrid$uby[-1]))[pgrid$p.in.grid,]
	sw <- as.matrix(expand.grid(pgrid$ubx[-1], pgrid$uby[-1]))[pgrid$p.in.grid,]
	
	if(!add)
	{
		plot(upx, upy, type = "n", ...)
	}
	
	if(length(set > 0))
	{
		for(i in set)
		{
			#Joins the x-coords of the border of the pixel
			x.border=c(nw[i,1], ne[i,1], se[i,1], sw[i,1], nw[i,1])
			#Joins the y-coords of the border of the pixel
			y.border=c(nw[i,2], ne[i,2], se[i,2], sw[i,2], nw[i,2])
			polygon(x.border,y.border, col=col, border=col)
		}
	}
}

#Check whether confidence set contains true clines
exceedance.success <- function(rejection.set, true.exceedance)
{
	overlap <- intersect(true.exceedance, rejection.set)
	success <- (length(overlap)==0)
	out <- NULL
	out$success <- success
	out$n.failed <- length(overlap)
	out$percent.contain=1-length(overlap)/length(true.exceedance)
	out
}

statistic.cv <- function(statistic.sim.obj, conf.level = .95)
{
	n <- length(statistic.sim.obj$statistic.sim)
	alternative <- statistic.sim.obj$alternative
	
	#Determine whether we can take the ith element of the sorted
	#statistics as the quantile, or if we need to average elements.  
	#No average is type == 1, average is type == 2.
	if((n * conf.level) == floor(n * conf.level))
	{
		type <- 1
	}else
	{
		type <- 2
	}
	
	#Determine critical value for statistic
	if(alternative == "less")
	{
		cv <- quantile(statistic.sim.obj$statistic.sim, 
			prob = 1 - conf.level, type = type)
	}else if(alternative == "greater")
	{
		cv <- quantile(statistic.sim.obj$statistic.sim, 
			prob = conf.level, type = type)
	}
	return(cv)
}

#Construct null or rejection region based on observed statistic
#and simulated statistic
exceedance.ci <- function(statistic.sim.obj, conf.level = .95, type = "null")
{
	alternative <- statistic.sim.obj$alternative
	cv <- statistic.cv(statistic.sim.obj, conf.level = conf.level)
	if(alternative == "less")
	{
		if(type == "null")
		{
			set <- which(statistic.sim.obj$statistic >= cv)
		}else
		{
			set <- which(statistic.sim.obj$statistic < cv)		
		}
	}else if(alternative == "greater")
	{
		if(type == "null")
		{
			set <- which(statistic.sim.obj$statistic <= cv)
		}else
		{
			set <- which(statistic.sim.obj$statistic > cv)		
		}
	}
	return(set)
}

#Graphically "activates" the pixels that a contour line intersects
#plot.pixel=function(which.pixel, partition, plot=FALSE, col="red")
#{
#  if(plot==TRUE) #If we need to plot the domain
#  {
#    xmin=min(partition$x.edges)  #Find the dimensions
#  	xmax=max(partition$x.edges)  #of the domain
#    ymin=min(partition$y.edges)
#    ymax=max(partition$y.edges)
#    
#    #Plot a blank domain
#    plot(c(xmin,xmax), c(ymin,ymax), type="n", xlab="x", ylab="y")
#  }
#  for(i in which.pixel)
#	{
#   #Constructs the y-coords of the border of the pixel
#	 x.border=c(partition$vert.coords[i,1], partition$vert.coords[i,2],
#              partition$vert.coords[i,2], partition$vert.coords[i,1],
#              partition$vert.coords[i,1])
#   #Constructs the x-coords of the border of the pixel
#	 y.border=c(partition$vert.coords[i,3], partition$vert.coords[i,3],
#              partition$vert.coords[i,4], partition$vert.coords[i,4],
#              partition$vert.coords[i,3])
#	 polygon(x.border,y.border,col=col,border=col)
#	}
#}
#
