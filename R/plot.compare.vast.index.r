

#' Plot comparison of multiple fit.vast model fits.
#' 
#' @param species String denoting the species of the analysis
#' @param vast.index.list Named list of the different models to compare. Each entry of the list should be an output of vast.frq.index.
#' @param which.idx Which index to plot either "idx.std" or "idx.frq"
#' @param region.names What should these regions be called? A character vector of names please
#' @param smooth.span Degree of smoothing used. Value between 0 and 1, where closer to 1 is more smooth.
#' @param pt.alpha alpha transparency of quarterly points or uncertainty band.
#' @param plot.uncertainty TRUE or FALSE: Plot the uncertainty. Not available for idx.frq type.
#' @param scale.type Either have yaxis fixed by panel ("fixed") or free ("free_y")
#' @param color.palette Color palette used to differentiate betweeen models. Give exact number of colors otherwise a colorRamp will be used.
#' @param save.dir Path to the directory where the outputs will be saved
#' @param save.name Name stem for the output, useful when saving many model outputs in the same directory
#' @export
#' @import magrittr
#' @importFrom data.table as.data.table
#' @importFrom data.table melt
#' @importFrom ggthemes theme_few
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 geom_hline

plot.compare.vast.index = function(species,vast.index.list,which.idx="idx.std",region.names = paste0("R",1:9),smooth.span=0.1,pt.alpha=0.1,plot.uncertainty=FALSE,scale.type = "fixed",color.palette=c("royalblue3","deepskyblue1","gold","orange1","indianred1","firebrick2","#AC2020"),save.dir,save.name)
{
	if(length(vast.index.list)<2)
	{
		stop("This is a comparison function so it needs more than 1 thing to compare, add another vast.frq.index object")
	}

	# merge into one data.table for plotting
	if(which.idx == "idx.std"){err.stem = "cv"}else{err.stem= "penwt"}
	vast.index.list = lapply(vast.index.list,function(x) x[[which(names(x)==which.idx)]])
	tmp.dt = data.table::as.data.table(vast.index.list[[1]]) %>% .[,model:=names(vast.index.list)[1]]
	tmp.dt.idx = tmp.dt[,c("yrqtr","model",region.names),with=FALSE] %>% data.table::melt(.,id.vars=c("yrqtr","model"),variable.name="region",value.name="index")
	tmp.dt.err = tmp.dt[,c("yrqtr","model",paste0(region.names,err.stem)),with=FALSE] %>% data.table::melt(.,id.vars=c("yrqtr","model"),variable.name="region",value.name="err") %>% .[,region:=gsub(err.stem,"",region)]
	plot.dt = merge(tmp.dt.idx,tmp.dt.err,by=c("yrqtr","model","region"))
	rm(list=c("tmp.dt","tmp.dt.idx","tmp.dt.err"))

	for(i in 2:length(vast.index.list))
	{
		tmp.dt = data.table::as.data.table(vast.index.list[[i]]) %>% .[,model:=names(vast.index.list)[i]]
		tmp.dt.idx = tmp.dt[,c("yrqtr","model",region.names),with=FALSE] %>% data.table::melt(.,id.vars=c("yrqtr","model"),variable.name="region",value.name="index")
		tmp.dt.err = tmp.dt[,c("yrqtr","model",paste0(region.names,err.stem)),with=FALSE] %>% data.table::melt(.,id.vars=c("yrqtr","model"),variable.name="region",value.name="err") %>% .[,region:=gsub(err.stem,"",region)]
		tmp.dt = merge(tmp.dt.idx,tmp.dt.err,by=c("yrqtr","model","region"))
		plot.dt = rbind(plot.dt,tmp.dt)
		rm(list=c("tmp.dt","tmp.dt.idx","tmp.dt.err"))
	}


		if(plot.uncertainty)
		{
			if(which.idx == "idx.std")
			{
				g = plot.dt %>% .[,model:=factor(as.character(model),levels=names(vast.index.list))] %>% .[,yr:=floor(yrqtr)] %>% .[,upper := exp(log(index)+2*err)] %>% .[,lower := exp(log(index)-2*err)] %>%
				.[,lapply(.SD,mean,na.rm=TRUE),by=.(model,region,yr),.SDcols=c("index","upper","lower")] %>%
				ggplot2::ggplot() + ggthemes::theme_few() + ggplot2::geom_hline(yintercept = 1,color="gray70") +
				ggplot2::xlab("Year") +
				ggplot2::ylab("CPUE") +
				ggplot2::ggtitle(paste0(species," standardized CPUE: vast indices (",which.idx,")")) +
				ggplot2::scale_color_manual(name="Model", values=colorRampPalette(color.palette)(length(vast.index.list))) +
				ggplot2::scale_fill_manual(name="Model", values=colorRampPalette(color.palette)(length(vast.index.list))) +
				ggplot2::geom_ribbon(ggplot2::aes(x=yr,ymin=lower,ymax=upper,fill=model),alpha=pt.alpha,color=NA) +
				ggplot2::geom_line(ggplot2::aes(x=yr,y=index,color=model),size=1.25) +
				ggplot2::facet_wrap(~region,drop=FALSE,scales = scale.type)
			} else {
				stop("Unable to compare uncertainty for the frq formatted index")
			}
		} else {
			g = plot.dt %>% .[,model:=factor(as.character(model),levels=names(vast.index.list))] %>%
			ggplot2::ggplot() + ggthemes::theme_few() +
			ggplot2::xlab("Year") +
			ggplot2::ylab("CPUE") +
			ggplot2::ggtitle(paste0(species," standardized CPUE: vast indices (",which.idx,")")) +
			ggplot2::scale_color_manual(name="Model", values=colorRampPalette(color.palette)(length(vast.index.list))) +
			ggplot2::geom_point(ggplot2::aes(x=yrqtr,y=index,color=model),alpha=pt.alpha) +
			ggplot2::geom_smooth(ggplot2::aes(x=yrqtr,y=index,color=model),se=FALSE,span=smooth.span) +
			ggplot2::facet_wrap(~region,drop=FALSE,scales = scale.type)
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