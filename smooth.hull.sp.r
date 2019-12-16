

smooth.hull.sp = function(pts,crs.ll,buffer.ll,vertices = 100, k=3)
# pts is a list
# inspired by https://gis.stackexchange.com/questions/24827/smoothing-polygons-in-contour-map/24929#24929
# https://stackoverflow.com/questions/42630703/create-buffer-around-spatial-data-in-r/42641283
# https://stackoverflow.com/questions/13577918/plotting-a-curve-around-a-set-of-points
{
	library(sp,quietly=TRUE)
	library(rgeos,quietly=TRUE)

	xy = as.matrix(as.data.frame(lapply(pts,"[",chull(pts)))) 
 	# Wrap k vertices around each end.
    n <- dim(xy)[1]
    if (k >= 1) {
        data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
    } else {
        data <- xy
    }

    # Spline the x and y coordinates.
    data.spline <- spline(1:(n+2*k), data[,1], n=vertices)
    x <- data.spline$x
    x1 <- data.spline$y
    x2 <- spline(1:(n+2*k), data[,2], n=vertices)$y

    # Retain only the middle part.
    smooth = cbind(x1, x2)[k < x & x <= n+k, ]
    smooth = rbind(smooth,smooth[1,])

	p = Polygon(smooth)
	ps = Polygons(list(p),1)
	sps = SpatialPolygons(list(ps))
	proj4string(sps) = crs.ll
	buff = gBuffer(sps, width=buffer.ll, quadsegs=10)
	return(buff)
}