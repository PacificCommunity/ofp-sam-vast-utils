

#' Plot DHARMa style residuals for VAST model.
#' 
#' @param vast.output The output from a fit.vast function call
#' @param n.samp Integer giving the number of records to plot the QQ relationship. Default is to use everyting.
#' @param seed The random seed for the subsampling
#' @param plot.type The type of metric to plot 'all', 'time', or 'knot', if you want to plot QQ by factor.
#' \describe{
#'   \item{'all'}{Observed nominal cpue by knot and timeblock}
#'   \item{'time'}{Predicted density cpue by knot and timeblock}
#'   \item{'knot'}{Predicted encounter rate by knot and timeblock}
#' }
#' @param palette.cols A color palette to use for plotting
#' @param smooth.span Value between 0 and 1. Degree of smoothing to use for the QQ plot
#' @param ts.start Integer giving the first year of the model
#' @param ts.step The timestep of the model 1 for year and 0.25 for quarter
#' @param save.dir Path to the directory where the outputs will be saved
#' @param save.name Name stem for the output, useful when saving many model outputs in the same directory
#' @export
#' @import magrittr
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom ggthemes theme_few
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_abline
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_colour_gradientn
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_sf
#' @importFrom ggplot2 coord_sf
#' @importFrom sf st_as_sf



plot.dharma.resid = function(vast.output,n.samp=999999,seed=123,plot.type="all",palette.cols = c("royalblue3","deepskyblue1","gold","orange1","indianred1","firebrick2","#AC2020"),smooth.span=0.1,ts.start=1952,ts.step=0.25,factor.grouping=5,save.dir,save.name)
{
	set.seed(seed)
	d.all = vast.output$dharma.all
	t_i = as.vector(vast.output$TmbData$t_iz) + 1 # ts vec
	k_i = data.table::data.table(knot=vast.output$Spatial_List$knot_i,lon=vast.output$Spatial_List$latlon_i[,2],lat=vast.output$Spatial_List$latlon_i[,1]) # knot vec
	n.obs = d.all$nObs
	if (! dir.exists(save.dir))dir.create(save.dir,recursive=TRUE)


	if(plot.type == "all")
	{
		if(n.samp<n.obs)
		{
			which.keep = sample(1:n.obs,n.samp)
			d.all$nObs = n.samp
			d.all$simulatedResponse = d.all$simulatedResponse[which.keep,]
			d.all$observedResponse = d.all$observedResponse[which.keep]
			d.all$scaledResiduals = d.all$scaledResiduals[which.keep]
			d.all$fittedPredictedResponse = d.all$fittedPredictedResponse[which.keep]
			t_i = t_i[which.keep]
			k_i = k_i[which.keep,]
		}

		g.dt = data.table::data.table(expected=sort((1:d.all$nObs)/(d.all$nObs + 1)),observed=sort(d.all$scaledResiduals))
		tmp.ks = suppressWarnings(ks.test(g.dt$observed, 'punif', alternative = "two.sided"))
		test.dt = data.table::data.table(p.value=round(tmp.ks$p.value, digits = 3),deviation=ifelse(tmp.ks$p.value < 0.05, "significant", "n.s."))
		test.dt$label = paste0("KS p-value: ",test.dt$p.value)

		g =g.dt %>% 
		ggplot2::ggplot() + ggthemes::theme_few() +
		ggplot2::xlab("Expected") +
		ggplot2::ylab("Observed") +
		ggplot2::ggtitle(paste0("QQ plot: ",round((d.all$nObs/n.obs)*100,digits=2),"% of records")) +
		ggplot2::geom_abline(slope=1, intercept=0, color="gray70",size=1.15) +
		ggplot2::geom_point(ggplot2::aes(x=expected,y=observed),color="blue") + 
		# ggplot2::geom_smooth(ggplot2::aes(x=expected,y=observed),se=FALSE,span=smooth.span,size=0.25) + 
		ggplot2::geom_text(size = 3,data = test.dt,mapping = ggplot2::aes(x = 0, y = 1, label = label),hjust=0,vjust=1)
		ggplot2::ggsave(paste0(save.name,"-QQunif.png"),plot=g, device = "png", path = save.dir,scale = 1, width = 16, height = 9, units = c("in"))


	} else if(plot.type == "time"){
		dt = data.table::data.table(factor = t_i,scaled.residual=d.all$scaledResiduals)
		u.factor = sort(unique(dt$factor))
		n.factor = length(u.factor)
		tmp.list = as.list(rep(NA,n.factor))
		for(i in 1:n.factor)
		{
			tmp.dt = dt[factor == u.factor[i]]
			tmp.n = nrow(tmp.dt)
			tmp.list[[i]] = data.table::data.table(factor=rep(u.factor[i],tmp.n),expected=sort((1:tmp.n)/(tmp.n + 1)),observed=sort(tmp.dt$scaled.residual))
			rm(list=c("tmp.dt","tmp.n"))
		}
		
		ts.vec = seq(from=ts.start,by=ts.step,length.out=n.factor)
		g.dt = data.table::rbindlist(tmp.list) %>% .[,factor:=ts.vec[factor]] %>% .[,factor.panel:=floor(factor/factor.grouping)*factor.grouping]

		u.factor.groups = sort(unique(g.dt$factor.panel))
		n.factor.groups = length(u.factor.groups)
		test.list = as.list(rep(NA,n.factor.groups))
		for(i in 1:n.factor.groups)
		{
			tmp.dt = g.dt[factor.panel == u.factor.groups[i]]
			tmp.ks = suppressWarnings(ks.test(tmp.dt$observed, 'punif', alternative = "two.sided"))
			test.list[[i]] = data.table::data.table(factor.panel=u.factor.groups[i],p.value=round(tmp.ks$p.value, digits = 3),deviation=ifelse(tmp.ks$p.value < 0.05, "significant", "n.s."))
			test.list[[i]]$label = paste0("KS p-value: ",test.list[[i]]$p.value)
			rm(list=c("tmp.dt","tmp.ks"))
		}
		test.dt = data.table::rbindlist(test.list)

		g =g.dt %>% 
		ggplot2::ggplot() + ggthemes::theme_few() +
		ggplot2::xlab("Expected") +
		ggplot2::ylab("Observed") +
		ggplot2::ggtitle("QQ plot by model timestep") +
		ggplot2::geom_abline(slope=1, intercept=0, color="gray70",size=1.15) + ggplot2::facet_wrap(~factor.panel) +
		# ggplot2::geom_point(ggplot2::aes(x=expected,y=observed,color=factor)) + 
		ggplot2::geom_smooth(ggplot2::aes(x=expected,y=observed,group=factor,color=factor),se=FALSE,span=smooth.span,size=0.25) + 
		ggplot2::scale_colour_gradientn(paste0("Factor type = ",plot.type),colors=palette.cols) + 
		ggplot2::geom_text(size = 3,data = test.dt,mapping = ggplot2::aes(x = 0, y = 1, label = label),hjust=0,vjust=1)
		ggplot2::ggsave(paste0(save.name,"-QQunif.time.png"),plot=g, device = "png", path = save.dir,scale = 1, width = 16, height = 9, units = c("in"))

		test.list = as.list(rep(NA,n.factor))
		for(i in 1:n.factor)
		{
			tmp.dt = g.dt[factor == ts.vec[u.factor[i]]]
			tmp.ks = suppressWarnings(ks.test(tmp.dt$observed, 'punif', alternative = "two.sided"))
			test.list[[i]] = data.table::data.table(factor.panel=ts.vec[u.factor[i]],p.value=round(tmp.ks$p.value, digits = 3),deviation=ifelse(tmp.ks$p.value < 0.05, "significant", "n.s."))

			rm(list=c("tmp.dt","tmp.ks"))
		}
		test.factor.dt = data.table::rbindlist(test.list) %>% .[,deviation:=factor(deviation,levels=c("n.s.","significant"))]
		g2 =test.factor.dt %>% 
		ggplot2::ggplot() + ggthemes::theme_few() +
		ggplot2::xlab("Year") +
		ggplot2::ylab("Observed") +
		ggplot2::ggtitle("QQ plot by model timestep") +
		ggplot2::geom_hline(yintercept=0.05, color="gray70",size=1.15) +
		ggplot2::geom_point(ggplot2::aes(x=factor.panel,y=p.value,color=deviation)) +
		ggplot2::scale_color_manual("Deviation from uniform",values=c("black","red"))
		ggplot2::ggsave(paste0(save.name,"-QQdev.time.png"),plot=g2, device = "png", path = save.dir,scale = 1, width = 16, height = 9, units = c("in"))


	} else if(plot.type == "knot"){
		dt = k_i
		dt$scaled.residual = d.all$scaledResiduals
		colnames(dt)[1] = "factor"
		u.factor = sort(unique(dt$factor))
		n.factor = length(u.factor)
		tmp.list = as.list(rep(NA,n.factor))
		for(i in 1:n.factor)
		{
			tmp.dt = dt[factor == u.factor[i]][order(scaled.residual)]
			tmp.n = nrow(tmp.dt)
			colnames(tmp.dt)[4] = "observed"
			tmp.dt$expected = sort((1:tmp.n)/(tmp.n + 1))
			tmp.list[[i]] = tmp.dt
			rm(list=c("tmp.dt","tmp.n"))
		}
		g.dt = data.table::rbindlist(tmp.list)
		n.factor.groups = round(n.factor/factor.grouping)
		factor.group.dt = g.dt[,lapply(.SD,mean),by=factor,.SDcols=c("lon","lat")]
		km = kmeans(as.data.frame(factor.group.dt[,.(lon,lat)]), n.factor.groups, iter.max = 100, nstart = 30)
		km.dt = data.table::data.table(round(km$centers)) %>% .[,old.cluster := 1:nrow(km$centers)] %>% .[order(lat,lon)] %>% .[,new.cluster := 1:nrow(km$centers)] %>% .[order(old.cluster)] %>% .[,factor.panel:=paste0("(",lon,", ",lat,")")]
		factor.group.dt$factor.panel = km.dt$factor.panel[km$cluster]
		g.dt = merge(g.dt[,.(factor,observed,expected)],factor.group.dt[,.(factor,factor.panel)],by="factor")

		u.factor.groups = sort(unique(g.dt$factor.panel))
		test.list = as.list(rep(NA,n.factor.groups))
		for(i in 1:n.factor.groups)
		{
			tmp.dt = g.dt[factor.panel == u.factor.groups[i]]
			tmp.ks = suppressWarnings(ks.test(tmp.dt$observed, 'punif', alternative = "two.sided"))
			test.list[[i]] = data.table::data.table(factor.panel=u.factor.groups[i],p.value=round(tmp.ks$p.value, digits = 3),deviation=ifelse(tmp.ks$p.value < 0.05, "significant", "n.s."))
			test.list[[i]]$label = paste0("KS p-value: ",test.list[[i]]$p.value)
			rm(list=c("tmp.dt","tmp.ks"))
		}
		test.dt = data.table::rbindlist(test.list)

		g = g.dt %>% 
		ggplot2::ggplot() + ggthemes::theme_few() +
		ggplot2::xlab("Expected") +
		ggplot2::ylab("Observed") +
		ggplot2::ggtitle("QQ plot by knots") +
		ggplot2::geom_abline(slope=1, intercept=0, color="gray70",size=1.15) + ggplot2::facet_wrap(~factor.panel) +
		# ggplot2::geom_point(ggplot2::aes(x=expected,y=observed,color=factor)) + 
		ggplot2::geom_smooth(ggplot2::aes(x=expected,y=observed,group=factor,color=factor),se=FALSE,span=smooth.span,size=0.25) + 
		ggplot2::scale_colour_gradientn(paste0("Factor type = ",plot.type),colors=palette.cols) + 
		ggplot2::geom_text(size = 3,data = test.dt,mapping = ggplot2::aes(x = 0, y = 1, label = label),hjust=0,vjust=1)
		ggplot2::ggsave(paste0(save.name,"-QQunif.space.png"),plot=g, device = "png", path = save.dir,scale = 1, width = 16, height = 9, units = c("in"))


		test.list = as.list(rep(NA,n.factor))
		for(i in 1:n.factor)
		{
			tmp.dt = g.dt[factor == u.factor[i]]
			tmp.ks = suppressWarnings(ks.test(tmp.dt$observed, 'punif', alternative = "two.sided"))
			test.list[[i]] = data.table::data.table(factor.panel=u.factor[i],p.value=round(tmp.ks$p.value, digits = 3),deviation=ifelse(tmp.ks$p.value < 0.05, "significant", "n.s."))
			test.list[[i]]$lon=vast.output$Spatial_List$latlon_x[u.factor[i],2]
			test.list[[i]]$lat=vast.output$Spatial_List$latlon_x[u.factor[i],1]
			rm(list=c("tmp.dt","tmp.ks"))
		}
		test.factor.dt = data.table::rbindlist(test.list) %>% .[lon<0,lon:=lon+360] %>% .[,deviation:=factor(deviation,levels=c("n.s.","significant"))]

		data("coast")
		data("wcpo.10N.shp")
		g2 = test.factor.dt %>% 
 		ggplot2::ggplot() + ggthemes::theme_few() +
		ggplot2::xlab("Longitude") +
		ggplot2::ylab("Latitude") +
		ggplot2::ggtitle("QQ residual uniformity by knot") +
		ggplot2::geom_point(ggplot2::aes(x=lon,y=lat,color=deviation),size=2) + 
		ggplot2::scale_color_manual("Deviation from uniform",values=c("black","red")) +
		ggplot2::geom_sf(data=sf::st_as_sf(coast)) + ggplot2::geom_sf(data=sf::st_as_sf(wcpo.10N.shp),fill=NA) + ggplot2::coord_sf(xlim = c(95, 295), ylim = c(-50, 60), expand = FALSE) 
 		ggplot2::ggsave(paste0(save.name,"-QQdev.space.png"),plot=g2, device = "png", path = save.dir,scale = 1, width = 16, height = 9, units = c("in"))


	} else {
		stop("Bad plot type. Try with 'all', 'time', or 'knot'")
	}

}

