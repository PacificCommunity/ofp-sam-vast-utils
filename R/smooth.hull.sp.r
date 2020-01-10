
#' Calculate a smooth hull around a set of points
#' Calculate a smooth polygon around a set of points using Bezier curves, with the smooth polygon passing through the vertices of the convex hull containing the set of points
#' once bezier control points are added the polygon is smoothed using chaikin's corner cutting algorithm
#' add buffer if desired
#' inspired by the following answers from stackexchange
#' Dan H's answer: https://gis.stackexchange.com/questions/24827/smoothing-polygons-in-contour-map/24929#24929
#' John Hughes's answer: https://math.stackexchange.com/questions/656500/given-a-point-slope-and-a-distance-along-that-slope-easily-find-a-second-p
#' https://stackoverflow.com/questions/42630703/create-buffer-around-spatial-data-in-r/42641283
#' 
#' @param pts List/Data Frame with two entries: one for x and one for y
#' @param crs.ll Character string of the crs for the current lat-lon projections
#' @param buffer.ll Size of buffer in degrees of lat/lon
#' @param d ratio bounded by (0,0.5) and controls the tightness of smooth to the original chull points
#' @return A spatial polygon object
#' @importFrom sp Polygon
#' @importFrom sp Polygons
#' @importFrom sp SpatialPolygons
#' @importFrom sp proj4string
#' @importFrom rgeos gBuffer
#' @export

smooth.hull.sp = function(pts,crs.ll,buffer.ll,d=0.15)
# pts is a list
{
    # calculate chull
        xy = as.matrix(as.data.frame(lapply(pts,"[",chull(pts))))

    # define vector to contain original points and bezier control points
        x=xy[,1]
        y=xy[,2]
        points = as.data.frame(cbind(x=x,y=y))
        points = rbind(points,points[1,])
        points.new = as.data.frame(cbind(x=x,y=y))[1,]

    # define storage structure in polar coordinates
    # we need this to ensure that we add the next appropriate Bezier control point in clockwise order
        origin.new = as.data.frame(matrix(colMeans(unique(points)),nrow=1,ncol=2))
        colnames(origin.new) = c("x","y")

        # shift first point relative to the new origin
            ref.vector = points[1,] + origin.new

    # calculate 1st Bezier control point for first vertex 
        m = (points$y[2] - tail(points$y,n=2)[1])/(points$x[2] - tail(points$x,n=2)[1])
        r = sqrt(1+m^2)
        d = d.scalar*sqrt((tail(points$x,n=2)[1]-points$x[1])^2+(tail(points$y,n=2)[1]-points$y[1])^2)    

        if(is.finite(m))
        {
            bez.a = data.frame(x=points$x[1]+d/r,y=points$y[1]+d*m/r)
            bez.b = data.frame(x=points$x[1]-d/r,y=points$y[1]-d*m/r)
        } else {
            # define case where slope is infinite
            bez.a = data.frame(x=points$x[1],y=points$y[1]+d)
            bez.b = data.frame(x=points$x[1],y=points$y[1]-d)
        }
            bez.points = rbind(bez.a,bez.b)

            # select as the 1st Bezier control point as the point with the greatest angle from the initial point
                bez.angle = apply(bez.points,1,function(x)calc.clockwise.angle(as.vector(as.matrix(points[1,] - origin.new)), as.vector(as.matrix(x - origin.new))))
                bez.points = bez.points[which(bez.angle == max(bez.angle)),]

            # hold aside to append to the end of points.new
                last.new = bez.points
                rm(list=c("m","r","d","bez.a","bez.b","bez.points","bez.angle"))

    # calculate 2nd Bezier control point for first vertex and append to x.new and y.new
        m = (points$y[2] - tail(points$y,n=2)[1])/(points$x[2] - tail(points$x,n=2)[1])
        r = sqrt(1+m^2)
        d = d.scalar*sqrt((points$x[2]-points$x[1])^2+(points$y[2]-points$y[1])^2)    

        if(is.finite(m))
        {
            bez.a = data.frame(x=points$x[1]+d/r,y=points$y[1]+d*m/r)
            bez.b = data.frame(x=points$x[1]-d/r,y=points$y[1]-d*m/r)
        } else {
            # define case where slope is infinite
            bez.a = data.frame(x=points$x[1],y=points$y[1]+d)
            bez.b = data.frame(x=points$x[1],y=points$y[1]-d)
        }
            bez.points = rbind(bez.a,bez.b)

            # select as the 2nd Bezier control point as the point with the smallest angle from the first point
                bez.angle = apply(bez.points,1,function(x)calc.clockwise.angle(as.vector(as.matrix(points[1,] - origin.new)), as.vector(as.matrix(x - origin.new))))
                bez.points = bez.points[which(bez.angle == min(bez.angle)),]

            # append to x.new & y.new
                points.new = rbind(points.new,bez.points)
                rm(list=c("m","r","d","bez.a","bez.b","bez.points","bez.angle"))

    # iterate across all remaining vertices & calculate Bezier control points for each vertex
        for(i in 2:(nrow(points)-1))
        {
            # define target vertex (V) and adjacent vertices
                V.a = points[i-1,]
                V = points[i,]
                V.b = points[i+1,]

            # define m, r, da & db
                m = (V.b$y - V.a$y)/(V.b$x - V.a$x)
                r = sqrt(1+m^2)
                da = d.scalar*sqrt((V.a$x-V$x)^2+(V.a$y-V$y)^2)    
                db = d.scalar*sqrt((V.b$x-V$x)^2+(V.b$y-V$y)^2)

            # 1st Bezier control point
                if(is.finite(m))
                {
                    bez.a = data.frame(x=V$x+da/r,y=V$y+da*m/r)
                    bez.b = data.frame(x=V$x-da/r,y=V$y-da*m/r)
                } else {
                    # define case where slope is infinite
                    bez.a = data.frame(x=V$x,y=V$y+da)
                    bez.b = data.frame(x=V$x,y=V$y-da)
                }
                    bez.points = rbind(bez.a,bez.b)
                # select as the 1st Bezier control point as the point closer to V.a
                    bez.angle = apply(bez.points,1,function(x)calc.clockwise.angle(as.vector(as.matrix(points[1,] - origin.new)), as.vector(as.matrix(x - origin.new))))
                    bez.points = bez.points[which(bez.angle == min(bez.angle)),]
                # append    
                    points.new = rbind(points.new,bez.points,V)
                    rm(list=c("bez.a","bez.b","bez.points","bez.angle"))

            # 2nd Bezier control point
                if(is.finite(m))
                {
                    bez.a = data.frame(x=V$x+db/r,y=V$y+db*m/r)
                    bez.b = data.frame(x=V$x-db/r,y=V$y-db*m/r)
                } else {
                    # define case where slope is infinite
                    bez.a = data.frame(x=V$x,y=V$y+db)
                    bez.b = data.frame(x=V$x,y=V$y-db)
                }
                    bez.points = rbind(bez.a,bez.b)
                # select as the 2nd Bezier control point as the point closer to V.b
                    bez.angle = apply(bez.points,1,function(x)calc.clockwise.angle(as.vector(as.matrix(points[1,] - origin.new)), as.vector(as.matrix(x - origin.new))))
                    bez.points = bez.points[which(bez.angle == max(bez.angle)),]
                # append    
                    points.new = rbind(points.new,bez.points)
                    rm(list=c("m","r","da","db","bez.a","bez.b","bez.points","bez.angle"))


        }
    # append last.new and first point 
        points.new = rbind(points.new,last.new,points[1,])
        rownames(points.new) = 1:nrow(points.new)

    # fit spline to each chunk of 4 points and append results (overlap vertices)
        start.sp = seq(from=1,to=nrow(points.new),by=3)[-length(seq(from=1,to=nrow(points.new),by=3))]
        end.sp = seq(from=1,to=nrow(points.new),by=3)[-1]

        smooth.pts = data.frame(x=NA,y=NA)
        for(i in 1:length(start.sp))
        {
            tmp.pts = as.matrix(points.new[start.sp[i]:end.sp[i],])
            for(j in 1:5)
            {
                tmp.pts = chaikin.open(tmp.pts)
            }
            colnames(tmp.pts) = c("x","y")
            tmp.pts = rbind(points.new[start.sp[i],],tmp.pts,points.new[end.sp[i],])
            smooth.pts = rbind(smooth.pts,tmp.pts)
            rm(list=c("tmp.pts"))

        }
        smooth.pts = na.omit(unique(smooth.pts))
        smooth.pts = rbind(smooth.pts,smooth.pts[1,])

        # convert to spatial polygon and buffer
            p = sp::Polygon(smooth.pts)
            ps = sp::Polygons(list(p),1)
            sps = sp::SpatialPolygons(list(ps))
            sp::proj4string(sps) = crs.ll
            buff = rgeos::gBuffer(sps, width=buffer.ll, quadsegs=10)
            return(buff)

        # plot(points,xlim=c(-3,3),ylim=c(-3,3),pch=16,cex=1.25)
        # points(origin.new,pch=16,col="hotpink")
        # lines(points,lty=3)
        # points(last.new,pch=16,col="green")
        # points(points.new,pch=16,col="red")
        # # lines(points.new,col="red")
        # lines(smooth.pts,col="green")


}
