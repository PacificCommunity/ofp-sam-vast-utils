
#' Function to convert points stored as eastings and northings to lat-lon
#' 
#' @param E Vector of eastings
#' @param N Vector of northings
#' @param crs.en Character string of the crs for the E-N projection
#' @param crs.ll Character string of the crs for the current lat-lon projections
#' @return Matrix of two columns: Lon and Lat
#' @importFrom sp coordinates
#' @importFrom sp proj4string
#' @importFrom sp CRS
#' @importFrom sp spTransform

#' @export

Convert_EN_to_LL_Fn.ndd = function(E, N, crs.en = "+proj=tpeqd +lat_1=0 +lon_1=155 +lat_2=0 +lon_2=209 +datum=WGS84 +ellps=WGS84 +units=km +no_defs",crs.ll = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# modified from Jim Thorson's FishStatUtils::Convert_LL_to_EastNorth_Fn.R
# pass x and y locations as Lon and Lat
# crs default is two point equidistant # OPTIONAL: +a=6005 appears to be the distance in km between the two points
{
  # Transform
	  dstart=data.frame(E=E, N=N) # that's the object
	  sp::coordinates(dstart) = c("E", "N")
	  sp::proj4string(dstart) = sp::CRS(crs.en) # that's the lat long projection
	  CRS.new = sp::CRS(crs.ll) # that's the eastings and northings projection
	  dstart.t = sp::spTransform(dstart, CRS.new) # here's where you transform

  # Clean up
	  dstart.t = cbind( "Lon"=dstart.t@coords[,"E"], "Lat"=dstart.t@coords[,'N'])
	  dstart.t[,1] = ifelse(dstart.t[,1]<0,dstart.t[,1]+360,dstart.t)
	  attr(dstart.t,"zone") = crs.ll

  # Return results
  	return( dstart.t )
}