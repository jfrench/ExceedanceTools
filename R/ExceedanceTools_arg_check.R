statistics_sim_arg_check <- function(obj, level, alternative)
{
	if(!is.null(obj))
	{
		if(is.null(obj$sim))
		{
			stop("krige.obj$sim cannot be NULL.  Try setting nsim > 0")
		}
	}
	if(!is.numeric(level) || length(level) > 1)
	{
		stop("level must be a numeric vector of length 1")
	}
	if(!(alternative == "less" || alternative == "greater" || alternative == "two.sided"))
	{
		stop('alternative must equal "two.sided" or "less" or "greater"')
	}
}