

#' Plot vast output by spatial knots and time block.
#' 
#' @param vast.output The output from a fit.vast function call
#' @param coast.shp supply a coast shapefile, if missing this defaults to the coast shapefile stored in the package
#' @param region.shp.list supply a shapefile with the regional structure, if missing this defaults to the WCPO 10N shapefile stored in the package
#' @param species String denoting the species of the analysis
#' @param plot.type The type of metric to plot 'pred.cpue', 'bin', 'pos', 'resid', 'bin.resid', 'pos.resid' or 'N'
#' \describe{
#'   \item{'obs.cpue'}{Observed nominal cpue by knot and timeblock}
#'   \item{'pred.cpue'}{Predicted density cpue by knot and timeblock}
#'   \item{'bin'}{Predicted encounter rate by knot and timeblock}
#'   \item{'pos'}{Predicted positive catch-rate by knot and timeblock}
#'   \item{'resid'}{Standardized residual by knot and timeblock}
#'   \item{'bin.resid'}{Standardized residual for encounter rate cpue by knot and timeblock}
#'   \item{'pos.resid'}{Standardized residual for positive catch-rate by knot and timeblock}
#'   \item{'N'}{Number of observations by knot and timeblock}
#' }
#' @param model.start.year An integer denoting the first year of data used in the model
#' @param t.block An integer denoting the number of years to aggregate together in time blocks
#' @param plot.obs.grid TRUE or FALSE. Plot the observations. Default is FALSE as otherwise the plot looks gross
#' @param x.extent Specify the xlim of the plot as c(xmin,xmax) in degrees 0 to 360
#' @param y.extent Specify the ylim of the plot as c(ymin,ymax) in degrees -90 to 90
#' @param scale.trans A character string denoting the type of ggplot2 scale transformation to use for the color scale. Built-in transformations include "asn", "atanh", "boxcox", "date", "exp", "hms", "identity", "log", "log10", "log1p", "log2", "logit", "modulus", "probability", "probit", "pseudo_log", "reciprocal", "reverse", "sqrt" and "time".
#' @param viridis.pal A character string indicating the colormap option to use. Four options are available: "magma" (or "A"), "inferno" (or "B"), "plasma" (or "C"), "viridis" (or "D", the default option) and "cividis" (or "E")
#' @param save.dir Path to the directory where the outputs will be saved
#' @param save.name Name stem for the output, useful when saving many model outputs in the same directory
#' @export
#' @import magrittr
#' @importFrom data.table as.data.table
#' @importFrom data.table data.table
#' @importFrom ggthemes theme_few
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_fill_viridis_c
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_sf
#' @importFrom ggplot2 coord_sf
#' @importFrom sf st_as_sf
#' @importFrom scales muted
#' @importFrom scales trans_new
#' @importFrom scales extended_breaks

plot.vast.space.time = function(vast.output,coast.shp,region.shp.list,species="BET",plot.type = "N",model.start.year=1952,t.block=4,plot.obs.grid=FALSE,x.extent=c(95, 295),y.extent=c(-50, 60),scale.trans="identity",virids.pal="D",save.dir,save.name)
{

	if(missing(region.shp.list))
	{
		data("wcpo.10N.shp")
		region.shp.list = wcpo.10N.shp
	}

	if(missing(coast.shp))
	{
		data("coast")
		coast.shp = coast
	}


	# pull out quantities and make objects
		Extrapolation_List = vast.output$Extrapolation_List
		Spatial_List = vast.output$Spatial_List
		TmbData = vast.output$TmbData
		Report = vast.output$Report
		D_xy = Report$D_gcy[,1,]
		R1_xy = Report$R1_gcy[,1,]
		R2_xy = Report$R2_gcy[,1,]

	# recreate Data_Geostat
		Data_Geostat = data.table::data.table(obs = TmbData$b_i, ts= as.vector(TmbData$t_iz)+1, knot = Spatial_List$knot_i, a_i = TmbData$a_i)
		colnames(Data_Geostat) = c("obs","ts","knot","a_i")
		Data_Geostat = cbind(Data_Geostat,Convert_EN_to_LL_Fn.ndd(Spatial_List$loc_i[,"E_km"], Spatial_List$loc_i[,"N_km"], crs.en = attr(Spatial_List$loc_i, "projargs"), crs.ll = attr(Spatial_List$loc_i, "origargs")))

	# add columns for yr.qtr, time.block
		yr.qtr.seq = seq(from=model.start.year,length.out=max(Data_Geostat$ts),by=0.25)
		time.block.seq = as.numeric(unique(floor(yr.qtr.seq/t.block)*t.block))
		Data_Geostat$yr.qtr = yr.qtr.seq[Data_Geostat$ts]
		Data_Geostat$time.block = floor(Data_Geostat$yr.qtr/t.block)*t.block
		Data_Geostat$qtr = ((Data_Geostat$yr.qtr %% 1)*4)+1
		Data_Geostat$obs.bin = ifelse(Data_Geostat$obs>0,1,0)

	# make extrapgrid.dt
		extrapgrid.dt=data.table::as.data.table(Extrapolation_List$Data_Extrap[,1:3])
		extrapgrid.dt$knot=as.vector(Spatial_List$NN_Extrap$nn.idx)

	# create dt by plot type
		if(plot.type == "obs.cpue")
		{
			dg = Data_Geostat[,.(metric=mean(obs)),by=.(time.block,knot)][order(time.block,knot)]
			dg$metric = dg$metric/mean(dg$metric)
			plot.title = "Nominal CPUE"
			legend.title = "Relative CPUE"

		}else if(plot.type == "pred.cpue"){
			Data_Geostat$knot.pred = sapply(1:nrow(Data_Geostat),function(x)D_xy[Data_Geostat$knot[x],Data_Geostat$ts[x]])
			dg = Data_Geostat[,.(metric=mean(knot.pred*a_i)),by=.(time.block,knot)][order(time.block,knot)]
			dg$metric = dg$metric/mean(dg$metric)
			plot.title = "Predicted Density"
			legend.title = "Relative density"
		}else if(plot.type == "bin"){
			Data_Geostat$knot.bin = sapply(1:nrow(Data_Geostat),function(x)R1_xy[Data_Geostat$knot[x],Data_Geostat$ts[x]])
			dg = Data_Geostat[,.(metric=mean(knot.bin)*100),by=.(time.block,knot)][order(time.block,knot)]
			plot.title = "Predicted Encouter Rate"
			legend.title = "Encounter rate (%)"
		}else if(plot.type == "pos"){
			Data_Geostat$knot.pos = sapply(1:nrow(Data_Geostat),function(x)R2_xy[Data_Geostat$knot[x],Data_Geostat$ts[x]])
			dg = Data_Geostat[,.(metric=mean(knot.pos*a_i)),by=.(time.block,knot)][order(time.block,knot)]
			dg$metric = dg$metric/mean(dg$metric)
			plot.title = "Predicted positive catch-rate"
			legend.title = "Relative catch-rate"
		}else if(plot.type == "resid"){
			Data_Geostat$knot.pred = sapply(1:nrow(Data_Geostat),function(x)D_xy[Data_Geostat$knot[x],Data_Geostat$ts[x]])
			Data_Geostat$pred.density = Data_Geostat$knot.pred*Data_Geostat$a_i
			dg = Data_Geostat[,.(metric=mean(obs-pred.density)),by=.(time.block,knot)][order(time.block,knot)]
			dg$metric = dg$metric/sd(dg$metric,na.rm=TRUE)
			plot.title = "Residual - Density"
			legend.title = "Std. Residual (Obs - Pred)"
		}else if(plot.type == "bin.resid"){
			Data_Geostat$knot.bin = sapply(1:nrow(Data_Geostat),function(x)R1_xy[Data_Geostat$knot[x],Data_Geostat$ts[x]])
			dg = Data_Geostat[,.(metric=mean(obs.bin-knot.bin)),by=.(time.block,knot)][order(time.block,knot)]
			dg$metric = dg$metric/sd(dg$metric,na.rm=TRUE)
			plot.title = "Residual - Encounter rate"
			legend.title = "Std. Residual (Obs - Pred)"
		}else if(plot.type == "pos.resid"){
			Data_Geostat$knot.pos = sapply(1:nrow(Data_Geostat),function(x)R2_xy[Data_Geostat$knot[x],Data_Geostat$ts[x]])
			Data_Geostat$pos.density = Data_Geostat$knot.pos*Data_Geostat$a_i
			dg = Data_Geostat[obs.bin==1,.(metric=mean(obs-pos.density)),by=.(time.block,knot)][order(time.block,knot)]
			dg$metric = dg$metric/sd(dg$metric,na.rm=TRUE)
			plot.title = "Residual - positive catch-rate"
			legend.title = "Std. Residual (Obs - Pred)"
		}else if(plot.type == "N"){
			dg = Data_Geostat[,.(metric=.N),by=.(time.block,knot)][order(time.block,knot)]
			plot.title = "Observations per knot"
			legend.title = "Number of observations"
		} else {
			stop("Bad plot type. Please use either: 'obs.cpue', 'pred.cpue', 'bin', 'pos', 'resid', 'bin.resid', 'pos.resid' or 'N'")
		}

		# make plot!

		tmp.dt = merge(extrapgrid.dt,dg,by="knot",allow.cartesian=TRUE)

		if(plot.obs.grid)
		{
			obs.points = Data_Geostat[,.(metric=.N),by=.(time.block,Lon,Lat)][order(time.block),.(time.block,Lon,Lat)]
			if(plot.type %in% c("N","obs.cpue","pred.cpue","bin","pos"))
			{
				g= tmp.dt %>%
				ggplot2::ggplot() + ggthemes::theme_few() +
				ggplot2::xlab("Longitude") +
				ggplot2::ylab("Latitude") +
				ggplot2::ggtitle(paste0(species," ",plot.title)) +
				ggplot2::geom_tile(ggplot2::aes(x=Lon,y=Lat,fill=metric)) + 
				ggplot2::geom_point(data=obs.points,ggplot2::aes(Lon,Lat),shape=1,size=0.5,color="white") +
				ggplot2::facet_wrap(~time.block,drop=FALSE) + 
				ggplot2::scale_fill_viridis_c(legend.title,option=virids.pal,trans=scale.trans) +
				ggplot2::geom_sf(data=sf::st_as_sf(coast.shp)) + ggplot2::geom_sf(data=sf::st_as_sf(region.shp.list),fill=NA) + 
				ggplot2::coord_sf(xlim = x.extent, ylim = y.extent, expand = FALSE)
			} else {
				# trim tails function comes from
				# https://stackoverflow.com/questions/44628130/ggplot2-dealing-with-extremes-values-by-setting-a-continuous-color-scale
				trim_tails <- function(range = c(-Inf, Inf)) scales::trans_new("trim_tails", 
	                transform = function(x) {
	                  force(range)
	                  desired_breaks <- scales::extended_breaks(n = 7)(x[x >= range[1] & x <= range[2]])
	                  break_increment <- diff(desired_breaks)[1]
	                  x[x < range[1]] <- range[1] - break_increment
	                  x[x > range[2]] <- range[2] + break_increment
	                  x
	                },
	                inverse = function(x) x,

	                breaks = function(x) {
	                  force(range)
	                  scales::extended_breaks(n = 7)(x)
	                },
	                format = function(x) {
	                  force(range)
	                  x[1] <- paste("<", range[1])
	                  x[length(x)] <- paste(">", range[2])
	                  x
	                })

				g= tmp.dt %>%
				ggplot2::ggplot() + ggthemes::theme_few() +
				ggplot2::xlab("Longitude") +
				ggplot2::ylab("Latitude") +
				ggplot2::ggtitle(paste0(species," ",plot.title)) +
				ggplot2::geom_tile(ggplot2::aes(x=Lon,y=Lat,fill=metric)) +
				ggplot2::geom_point(data=obs.points,ggplot2::aes(Lon,Lat),shape=1,size=0.5) +
				ggplot2::facet_wrap(~time.block,drop=FALSE) + 
				ggplot2::scale_fill_gradient2(legend.title,low = scales::muted("blue"),mid = "white",high = scales::muted("red"),trans = trim_tails(range = c(-3,3))) +
				ggplot2::geom_sf(data=sf::st_as_sf(coast.shp)) + ggplot2::geom_sf(data=sf::st_as_sf(region.shp.list),fill=NA) + 
				ggplot2::coord_sf(xlim = x.extent, ylim = y.extent, expand = FALSE)
			}
		} else {
			if(plot.type %in% c("N","obs.cpue","pred.cpue","bin","pos"))
			{
				g= tmp.dt %>%
				ggplot2::ggplot() + ggthemes::theme_few() +
				ggplot2::xlab("Longitude") +
				ggplot2::ylab("Latitude") +
				ggplot2::ggtitle(paste0(species," ",plot.title)) +
				ggplot2::geom_tile(ggplot2::aes(x=Lon,y=Lat,fill=metric)) + 
				ggplot2::facet_wrap(~time.block,drop=FALSE) + 
				ggplot2::scale_fill_viridis_c(legend.title,option=virids.pal,trans=scale.trans) +
				ggplot2::geom_sf(data=sf::st_as_sf(coast.shp)) + ggplot2::geom_sf(data=sf::st_as_sf(region.shp.list),fill=NA) + 
				ggplot2::coord_sf(xlim = x.extent, ylim = y.extent, expand = FALSE)
			} else {
				# trim tails function comes from
				# https://stackoverflow.com/questions/44628130/ggplot2-dealing-with-extremes-values-by-setting-a-continuous-color-scale
				trim_tails <- function(range = c(-Inf, Inf)) scales::trans_new("trim_tails", 
	                transform = function(x) {
	                  force(range)
	                  desired_breaks <- scales::extended_breaks(n = 7)(x[x >= range[1] & x <= range[2]])
	                  break_increment <- diff(desired_breaks)[1]
	                  x[x < range[1]] <- range[1] - break_increment
	                  x[x > range[2]] <- range[2] + break_increment
	                  x
	                },
	                inverse = function(x) x,

	                breaks = function(x) {
	                  force(range)
	                  scales::extended_breaks(n = 7)(x)
	                },
	                format = function(x) {
	                  force(range)
	                  x[1] <- paste("<", range[1])
	                  x[length(x)] <- paste(">", range[2])
	                  x
	                })

				g= tmp.dt %>%
				ggplot2::ggplot() + ggthemes::theme_few() +
				ggplot2::xlab("Longitude") +
				ggplot2::ylab("Latitude") +
				ggplot2::ggtitle(paste0(species," ",plot.title)) +
				ggplot2::geom_tile(ggplot2::aes(x=Lon,y=Lat,fill=metric)) + 
				ggplot2::facet_wrap(~time.block,drop=FALSE) + 
				ggplot2::scale_fill_gradient2(legend.title,low = scales::muted("blue"),mid = "white",high = scales::muted("red"),trans = trim_tails(range = c(-3,3))) +
				ggplot2::geom_sf(data=sf::st_as_sf(coast.shp)) + ggplot2::geom_sf(data=sf::st_as_sf(region.shp.list),fill=NA) + 
				ggplot2::coord_sf(xlim = x.extent, ylim = y.extent, expand = FALSE)
			}
		}

		# write.out
		if(!missing(save.dir))
		{
			if(missing(save.name))
			{
				stop("How can you save the output if you haven't specified the directory? Please specify save.dir.")
			} else {
				if (! dir.exists(save.dir))dir.create(save.dir,recursive=TRUE)
				ggplot2::ggsave(paste0(save.name,".png"),plot=g, device = "png", path = save.dir,scale = 1, width = 16, height = 9, units = c("in"))
			}
		} 
			
		return(g)
}
