
#' Get the area of each longitude by latitude cell. Adjust for change in area with latitude.
#' 
#' @param lond Vector of longitudes (midpoint of the cell)
#' @param latd Vector of latitudes (midpoint of the cell)
#' @param cell.size Integer: size of the cell in degrees
#' @param crs.ll Character string of the crs for the current lat-lon projections
#' @param crs.en.type What type of projection for calculating the area? 'tpeqd', 'mollweide', or 'albers equal area'. Note that the default 'albers equal area' is hard coded for the central pacific
#' @return Either a vector of areas in km2
#' @importFrom sp SpatialPolygons
#' @importFrom sp CRS
#' @importFrom rgeos gArea
#' @importFrom sp spTransform
#' @export
#' 
area.cell = function(lond,latd,cell.size=rep(1,length(lond)),crs.ll = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",crs.en.type = "albers equal area")
{
	if(crs.en.type=="tpeqd")
	{
		crs.en = paste0("+proj=tpeqd +lat_1=",mean(latd,na.rm=TRUE)," +lon_1=",round(min(lond,na.rm=TRUE) + (1/3)*abs(diff(range(lond,na.rm=TRUE))))," +lat_2=",mean(latd,na.rm=TRUE)," +lon_2=",round(min(lond,na.rm=TRUE) + (2/3)*abs(diff(range(lond,na.rm=TRUE))))," +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
	} else if(crs.en.type=="mollweide"){
		crs.en = paste0("+proj=moll +lon_0=",mean(lond,na.rm=TRUE)," +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs")
	} else if(crs.en.type=="albers equal area"){
		# crs.en = paste0("+proj=aea +lat_1=",round(min(latd,na.rm=TRUE) + (1/4)*abs(diff(range(latd,na.rm=TRUE))))," +lat_2=",round(min(latd,na.rm=TRUE) + (3/4)*abs(diff(range(latd,na.rm=TRUE))))," +lat_0=",round(min(latd,na.rm=TRUE) + (2/4)*abs(diff(range(latd,na.rm=TRUE))))," +lon_0=",mean(lond,na.rm=TRUE)," +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
		crs.en = paste0("+proj=aea +lat_1=5 +lat_2=45 +lat_0=25 +lon_0=210 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")

	} else {
		stop("Incorrect crs type is available now. Please set crs.en.type to 'tpeqd'or 'mollweide' or 'albers equal area' ")
	}

	coords = data.frame(x=lond,y=latd)
	tmp.fn = function(id,coords,crs.ll)
	{	
		coords = coords[id,]
		coords = data.frame(x=coords[[1]]+c(0.5*cell.size,-0.5*cell.size,-0.5*cell.size,0.5*cell.size),y=coords[[2]]+c(-0.5*cell.size,-0.5*cell.size,0.5*cell.size,0.5*cell.size))
		poly = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(thicken.poly(rbind(coords, coords[1, ]),100))), ID = id)),proj4string=sp::CRS(crs.ll))
		return(poly)
	}

	# area is invariant within the same latitude 
	u.lat.df = data.frame(x=rep(mean(coords$x,na.rm=TRUE),length(unique(coords$y))),y=sort(unique(coords$y)))
		
	poly.1x1 = lapply(1:nrow(u.lat.df),tmp.fn,coords=u.lat.df,crs.ll=crs.ll)
	poly.1x1 = do.call("rbind",poly.1x1)
	poly.1x1.trans = sp::spTransform(poly.1x1, crs.en)


	u.lat.df$area.km2 = rgeos::gArea(poly.1x1.trans,byid=TRUE)
	coords$area.km2 = u.lat.df$area.km2[match(coords$y,u.lat.df$y)]
	
	return(coords$area.km2)
	
}