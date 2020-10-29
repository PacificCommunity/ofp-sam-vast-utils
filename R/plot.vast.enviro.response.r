

#' Plot estimated catchability effects. Currently assumes only one set of interactions and no splines on catchability.
#' 
#' @param vast.output The output from a fit.vast function call
#' @param enviro.formula The formula for the environmental covariates (abundance) that was passed to fit.vast as a part of the enviro argument
#' @param report.enviro Named list ('mean', 'sd', 'min', and 'max') where each element of the list is a vector of length covariates in enviro.formula containing the appropriate metric. This is assuming that the environmental covariate was scaled to mean 0 and sd 1 prior to passing to fit.vast.
#' @param save.dir Path to the directory where the outputs will be saved
#' @param save.name Name stem for the output, useful when saving many model outputs in the same directory
#' @export
#' @import magrittr
#' @import splines
#' @importFrom data.table as.data.table
#' @importFrom data.table melt
#' @importFrom ggthemes theme_few
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_color_viridis_c
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 geom_hline

plot.vast.enviro.response = function(vast.output,enviro.formula,report.enviro,save.dir,save.name)
{
	# identify factors
		enviro.factors = as.vector(sapply(strsplit(strsplit(enviro.formula,"~")[[1]][2],"[:+()]+"),trimws))
		if(length(which(enviro.factors%in%c("bs","")))>0)
		{
			enviro.factors = enviro.factors[-which(enviro.factors%in%c("bs",""))]
		}

	# create dependent variable
		n.dep = 1000
		dep.scale = matrix(NA,nrow=n.dep,ncol=length(enviro.factors))
		dep = dep.scale
		for(j in 1:length(enviro.factors))
		{
			dep[,j] = seq(from=report.enviro$min[j],to=report.enviro$max[j],length.out=n.dep)
			dep.scale[,j] = (dep[,j] - report.enviro$mean[j])/report.enviro$sd[j]
		}
		colnames(dep.scale) = colnames(dep) = enviro.factors
		dep.scale = as.data.frame(dep.scale)
		dep = as.data.frame(dep)

	# recreate model.matrix
		X = model.matrix( update.formula(as.formula(enviro.formula), ~.+0), data=dep.scale )[,,drop=FALSE]

	# grab parameters
		Report = vast.output$Report
		par.1 = vast.output$Opt$par[grep("gamma1_ctp",names(vast.output$Opt$par))]
		par.2 = vast.output$Opt$par[grep("gamma2_ctp",names(vast.output$Opt$par))]
		names(par.2) = names(par.1) = colnames(X)

	# make response
		rep2 = rep1 = matrix(NA,nrow=n.dep,ncol=length(enviro.factors))
		colnames(rep2) =colnames(rep1) = enviro.factors

		for(j in 1:length(enviro.factors))
		{
			which.par = grep(enviro.factors[j],colnames(X))
			rep1[,j] = X[,which.par] %*% matrix(par.1[which.par],nrow=length(which.par),ncol=1)
			rep2[,j] = X[,which.par] %*% matrix(par.2[which.par],nrow=length(which.par),ncol=1)
		}
		rep1 = as.data.frame(rep1)
		rep2 = as.data.frame(rep2)

	# convert to data.table & combine for plotting
		dep = data.table::as.data.table(dep) %>% data.table::melt(.,measure.vars=enviro.factors,variable.name = "Covariate")
		dep = rbind(dep,dep)
		dep$Component = c(rep("Encounter probability",n.dep*length(enviro.factors)),rep("Positive catch",n.dep*length(enviro.factors)))

		rep1 = data.table::as.data.table(rep1) %>% data.table::melt(.,measure.vars=enviro.factors,variable.name = "Covariate",value.name="Response")
		rep1$Component = "Encounter probability"
		rep2 = data.table::as.data.table(rep2) %>% data.table::melt(.,measure.vars=enviro.factors,variable.name = "Covariate",value.name="Response")
		rep2$Component = "Positive catch"
		rep = rbind(rep1,rep2)

		g = dep %>% .[,Response:=rep$Response] %>% .[order(Component,Covariate,value)] %>%
			ggplot2::ggplot() + ggthemes::theme_few() + ggplot2::geom_hline(yintercept = 0,color="gray70") +
			ggplot2::xlab("Covariate") +
			ggplot2::ylab("Response") +
			ggplot2::ggtitle("Estimated response to environmental covariates") +
			ggplot2::geom_point(ggplot2::aes(x=value,y=Response,color=Response)) +
			ggplot2::scale_color_viridis_c("Estimated response",option="D") +
			ggplot2::facet_grid(Component~Covariate,drop=FALSE,scale="free_x")
			# write.out
		
		if(!missing(save.dir))
		{
			if(missing(save.name))
			{
				stop("How can you save the output if you haven't specified the directory? Please specify save.dir.")
			} else {
				if (! dir.exists(save.dir))dir.create(save.dir,recursive=TRUE)
				ggplot2::ggsave(paste0(save.name,".png"),plot=g, device = "png", path = save.dir,scale = 1, width = 9, height = 9, units = c("in"))
			}
		} 
		
		return(g)
}