
#' Function to make polygon "thicker", useful when transforming to a different crs
#' 
#' @param coords matrix of x and y coordinates
#' @param n number of points to "thicken" each segment by
#' @return Matrix of two columns: x and y
#' @export

thicken.poly = function(coords,n)
{
	thick.poly = c()
	for(i in 2:nrow(coords))
	{
		x.thick = seq(from=coords[i-1,1],to=coords[i,1],length.out=n)
		y.thick = seq(from=coords[i-1,2],to=coords[i,2],length.out=n)
		thick.poly = rbind(thick.poly,cbind(x.thick,y.thick))
	}
	return(thick.poly)
}


