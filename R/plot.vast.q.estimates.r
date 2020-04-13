

#' Plot estimated catchability effects. Currently assumes only one set of interactions and no splines on catchability.
#' 
#' @param vast.output The output from a fit.vast function call
#' @param coef.names Character vector of the factors used to create Q_ik
#' @param level.names Character vector of how to name the factors in the plot
#' @param error.structure Character string denoting the error structure used to fit the model. This is needed to backtransform the parameter estimates.
#' @param save.dir Path to the directory where the outputs will be saved
#' @param save.name Name stem for the output, useful when saving many model outputs in the same directory
#' @export
#' @import magrittr
#' @importFrom data.table as.data.table
#' @importFrom scales muted
#' @importFrom ggthemes theme_few
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 geom_hline
#' @importFrom boot inv.logit


plot.vast.q.estimates = function(vast.output, coef.names = c("flg.grp","syn.hbf.imb"),level.names = c("Flag group","HBF (pred)"),error.structure = "dln",save.dir,save.name)
# only for objects where Q_ik is estimated
# also assumes that there are no splines
# current only supports 1 interaction between factors!!
{
	TmbData = vast.output$TmbData
	Report = vast.output$Report
	Q_ik = TmbData$Q_ik

	if(error.structure != "dln")
	{
		stop("Not defined for models not using a delta-lognormal error structure")
	} else {
		par.1 = boot::inv.logit(vast.output$Opt$par[grep("lambda1_k",names(vast.output$Opt$par))])
		par.2 = exp(vast.output$Opt$par[grep("lambda2_k",names(vast.output$Opt$par))])
	}



		levels.list.1 = as.list(rep(NA,length(coef.names)))
		levels.list.2 = as.list(rep(NA,length(coef.names)))
		names(levels.list.2) = names(levels.list.1) = coef.names

		for(i in 1:length(coef.names))
		{
			names.tmp = colnames(Q_ik)[grep(coef.names[i],colnames(Q_ik))]
			lvls.tmp = sort(unique(unname(sapply(names.tmp,function(x)gsub(coef.names[i],"",strsplit(x,":")[[1]][grep(coef.names[i],strsplit(x,":")[[1]])],fixed=TRUE)))))
			par2.tmp = par1.tmp = rep(NA,length(lvls.tmp))
			names(par2.tmp) = names(par1.tmp) = lvls.tmp

			Q_ik2.tmp = Q_ik1.tmp = Q_ik
			for(j in 1:length(lvls.tmp))
			{
				cols.tmp = grep(paste0(coef.names[i],lvls.tmp[j]),names.tmp,fixed=TRUE)
				par1.tmp[j] = mean(par.1[cols.tmp],na.rm=TRUE)
				par2.tmp[j] = mean(par.2[cols.tmp],na.rm=TRUE)
				rm(list=c("cols.tmp"))
			}
			levels.list.1[[i]] = par1.tmp - mean(par1.tmp,na.rm=TRUE)
			levels.list.2[[i]] = par2.tmp - mean(par2.tmp,na.rm=TRUE)
			rm(list=c("names.tmp","lvls.tmp","par1.tmp","par2.tmp"))
		}


	
	# create data.table to store the different results
		dg.1 = data.table::as.data.table(data.frame(Component = rep("Encounter probability",length(levels.list.1[[1]])),Factor = rep(names(levels.list.1)[1],length(levels.list.1[[1]])),Level=names(levels.list.1[[1]]),Estimate=levels.list.1[[1]]))
		dg.2 = data.table::as.data.table(data.frame(Component = rep("Positive catch",length(levels.list.2[[1]])),Factor = rep(names(levels.list.2)[1],length(levels.list.2[[1]])),Level=names(levels.list.2[[1]]),Estimate=levels.list.2[[1]]))

		for(i in 2:length(coef.names))
		{
			dg.1 = rbind(dg.1,data.table::as.data.table(data.frame(Component = rep("Encounter probability",length(levels.list.1[[i]])),Factor = rep(names(levels.list.1)[i],length(levels.list.1[[i]])),Level=names(levels.list.1[[i]]),Estimate=levels.list.1[[i]])))
			dg.2 = rbind(dg.2,data.table::as.data.table(data.frame(Component = rep("Positive catch",length(levels.list.2[[i]])),Factor = rep(names(levels.list.2)[i],length(levels.list.2[[i]])),Level=names(levels.list.2[[i]]),Estimate=levels.list.2[[i]])))

		}

		g = rbind(dg.1,dg.2) %>% .[,Factor:=factor(as.character(Factor),levels = coef.names,labels=level.names)] %>%
			ggplot2::ggplot() + ggthemes::theme_few() + ggplot2::geom_hline(yintercept = 0,color="black") +
			ggplot2::xlab("Level") +
			ggplot2::ylab("Effect") +
			ggplot2::ggtitle("Relative estimated catchability effects") +
			ggplot2::scale_x_discrete(drop=TRUE) +
			ggplot2::geom_bar(ggplot2::aes(x=Level,y=Estimate,fill=Estimate),stat = "identity") +
			ggplot2::scale_fill_gradient2("Relative effect",low = scales::muted("blue"),mid = "gray90",high = scales::muted("red")) +
			ggplot2::facet_grid(Component~Factor,drop=TRUE, space="free_x",scale="free_x")
	
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

