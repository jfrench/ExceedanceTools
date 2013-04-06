create.pgrid <- function(xmin, xmax, ymin, ymax, nx, ny, midpoints = FALSE,
	poly.coords = NULL)
{
	if(midpoints)
	{
	
		xstep <- (xmax-xmin)/(nx - 1)	#Calculate the pixel width
		ystep <- (ymax-ymin)/(ny - 1)	#Calculate the pixel height
	

		#Determine x and y midpoints of all of the pixels	
		upx <- xmin + 0:nx * xstep 
		upy <- ymin + 0:nx * ystep 
		
		#Create boundaries for pixels
		ubx <- xmin + 0:(nx + 1) * xstep - xstep/2
		uby <- ymin + 0:(ny + 1) * ystep - ystep/2
			
	}
	else
	{
		xstep <- (xmax-xmin)/nx	#Calculate the pixel width
		ystep <- (ymax-ymin)/ny	#Calculate the pixel height

		#Determine x and y midpoints of all of the pixels	
		upx <- xmin + xstep/2 + 0:(nx-1) * xstep 
		upy <- ymin + ystep/2 + 0:(ny-1) * ystep 
		
		#Create boundaries for pixels
		ubx <- xmin + 0:nx * xstep
		uby <- ymin + 0:ny * ystep
	}

	
	#If coords are supplied, create pgrid based on whether points
	#are contained in polygon of poly.coords.
	if(!is.null(poly.coords))
	{
		all.grid <- as.matrix(expand.grid(upx, upy))
		#Determine points of rectangular grid (based on xgrid and ygrid)
		#within poly.coords
		pip <- inout(all.grid, poly.coords, bound = TRUE)
		#Extract prediction coordinates within border 
		pgrid <- as.matrix(all.grid[pip == 1,])
		#Determine number of prediction locations
		np <- nrow(pgrid)
		#Determine which points are a prediction coordinates within 
		#rectangular grid
		p.in.grid <- (pip == 1)
	}else
	{
		pgrid <- as.matrix(expand.grid(upx, upy))
		np <- length(upx) * length(upy)
		p.in.grid <- rep(TRUE, np)
	}
	
	out <- (list(pgrid = as.matrix(pgrid), np = np, p.in.grid = p.in.grid, 
		ubx = ubx, uby = uby, upx = upx, upy = upy))	
	class(out) <- "pgrid"
	return(out)
	
}

create.pgrid2 <- function(xgrid, ygrid, midpoints = FALSE, poly.coords = NULL)
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
	#are contained in polygon of poly.coords.
	if(!is.null(poly.coords))
	{
		all.grid <- as.matrix(expand.grid(upx, upy))
		#Determine points of rectangular grid (based on xgrid and ygrid)
		#within poly.coords
		pip <- inout(all.grid, poly.coords, bound = TRUE)
		#Extract prediction coordinates within border 
		pgrid <- as.matrix(all.grid[pip == 1,])
		#Determine number of prediction locations
		np <- nrow(pgrid)
		#Determine which points are a prediction coordinates within 
		#rectangular grid
		p.in.grid <- (pip == 1)
	}else
	{
		pgrid <- as.matrix(expand.grid(upx, upy))
		np <- length(upx) * length(upy)
		p.in.grid <- rep(TRUE, np)
	}
	
	out <- (list(pgrid = as.matrix(pgrid), np = np, p.in.grid = p.in.grid, 
		ubx = ubx, uby = uby, upx = upx, upy = upy))	
	class(out) <- "pgrid"
	return(out)
}

statistic.sim <- function(krige.obj, level, alternative = "less", ...)
{
	statistics_sim_arg_check(krige.obj, level, alternative)

	nsim <- ncol(krige.obj$sim)
	stat.sim <- numeric(nsim)
	if(alternative != "two.sided")
	{
		statistic <- (krige.obj$pred - level)/sqrt(krige.obj$mspe)
	}else
	{
		statistic <- abs(krige.obj$pred - level)/sqrt(krige.obj$mspe)
	}

	if(alternative == "less")
	{
		for(i in 1:nsim)
		{
			which.exceedance.sim <- which(krige.obj$sim[, i] >= level)
			if(length(which.exceedance.sim) > 0)
			{
				stat.sim[i] <- min(statistic[which.exceedance.sim])
			}
		}
	}else if(alternative == "greater")
	{
		for(i in 1:nsim)
		{
			which.exceedance.sim <- which(krige.obj$sim[, i] <= level)
			if(length(which.exceedance.sim) > 0)
			{
				stat.sim[i] <- max(statistic[which.exceedance.sim])
			}
		}
	}
	else
	{
		# Check arguments.  Create unspecified arguments if needed.
		arglist <- list(...)
		argnames <- names(arglist)
		user.cov <- arglist$user.cov
		pgrid <- arglist$pgrid
		npx <- length(pgrid$upx)
		npy <- length(pgrid$upy)

		for(m in 1:nsim)
		{
			simmat <- matrix(krige.obj$sim[, m], nrow = npx, ncol = npy)
			cL <- contourLines(pgrid$upx, pgrid$upy, simmat, levels = level)
			
			if(length(cL) > 0)
			{
				obj <- user.cov(cLcoords = get.contours(cL), ...)
				krige.cL <- krige.uk(arglist$y, obj$V, obj$Vp, obj$Vop, arglist$X, obj$Xp)
				stat.sim[m] <- max(abs(krige.cL$pred - level)/sqrt(krige.cL$mspe))
			}
		}
	}

	return(list(statistic = statistic, statistic.sim = stat.sim,  
		alternative = alternative, level = level))
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
	}else
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
	else
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

plot.pgrid <- function(x, set, col = "gray", add = FALSE, ...)
{
	if(is.null(x$ubx) || is.null(x$uby))
	{
		stop("pgrid must contain boundaries.  Use create.pgrid or create.pgrid2 function")
	}
	nx <- length(x$ubx)
	ny <- length(x$uby)
	
	#grid point bound for nw, ne, se, and sw of pcoords
	nw <- as.matrix(expand.grid(x$ubx[-1], x$uby[-ny]))[x$p.in.grid,]
	ne <- as.matrix(expand.grid(x$ubx[-nx], x$uby[-ny]))[x$p.in.grid,]
	se <- as.matrix(expand.grid(x$ubx[-nx], x$uby[-1]))[x$p.in.grid,]
	sw <- as.matrix(expand.grid(x$ubx[-1], x$uby[-1]))[x$p.in.grid,]
	
	if(!add)
	{
		plot(x$upx, x$upy, type = "n", ...)
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


