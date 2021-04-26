
#' Function to convert points stored as lat-lon to eastings and northings
#' 
#' @param Lon Vector of longitudes
#' @param Lat Vector of latitudes
#' @param crs.en Character string of the crs for the E-N projection
#' @param crs.ll Character string of the crs for the current lat-lon projections
#' @return Matrix of two columns: E_km and N_km
#' @importFrom sp coordinates
#' @importFrom sp proj4string
#' @importFrom sp CRS
#' @importFrom sp spTransform
#' @export

Convert_LL_to_EastNorth_Fn.ndd = function(Lon, Lat, crs.en = "+proj=tpeqd +lat_1=0 +lon_1=155 +lat_2=0 +lon_2=209 +datum=WGS84 +ellps=WGS84 +units=km +no_defs",crs.ll = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# modified from Jim Thorson's FishStatUtils::Convert_LL_to_EastNorth_Fn.R
# pass x and y locations as Lon and Lat
# crs default is two point equidistant # OPTIONAL: +a=6005 appears to be the distance in km between the two points
{

  # Transform
	  dstart=data.frame(lon=Lon, lat=Lat) # that's the object
	  sp::coordinates(dstart) = c("lon", "lat")
	  sp::proj4string(dstart) = sp::CRS(crs.ll) # that's the lat long projection
	  CRS.new = sp::CRS(crs.en) # that's the eastings and northings projection
	  dstart.t = sp::spTransform(dstart, CRS.new) # here's where you transform

  # Clean up
	  dstart.t = cbind( "E_km"=dstart.t@coords[,"lon"], "N_km"=dstart.t@coords[,'lat'])
	  attr(dstart.t,"zone") = crs.en

  # Return results
  	return( dstart.t )
}