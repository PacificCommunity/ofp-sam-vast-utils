
#' Influence plot for fitted model from fit.vast. Currently assumes only one set of interactions and no splines on catchability.
#' 
#' @param vast.output The output from a fit.vast function call
#' @param model.start.year An integer denoting the first year of data used in the model
#' @param coef.names Character vector of the factors used to create Q_ik
#' @param level.names Character vector of how to name the factors in the plot
#' @param error.structure Character string denoting the error structure used to fit the model. This is needed to backtransform the parameter estimates.
#' @param smooth.span Degree of smoothing used. Value between 0 and 1, where closer to 1 is more smooth.
#' @param pt.alpha alpha transparency of quarterly points or uncertainty band.
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
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 geom_hline
#' @importFrom boot inv.logit

plot.vast.influ = function(vast.output,model.start.year=1952, coef.names = c("flg.grp","syn.hbf.imb"),level.names = c("Flag group","HBF (pred)"),error.structure = "dln", pt.alpha = 0.2,smooth.span = 0.1,save.dir,save.name)
# only for objects where Q_ik is estimated
# current only supports 1 interaction between factors!!
{
	TmbData = vast.output$TmbData
	Report = vast.output$Report
	Q_ik = TmbData$Q_ik

	if(!(error.structure %in% c("dln","dg","dp","dnb")))
	{
		stop("Not defined for models not using logit and log link functions")
	} else {
		par.1 = boot::inv.logit(vast.output$Opt$par[grep("lambda1_k",names(vast.output$Opt$par))])
		par.2 = exp(vast.output$Opt$par[grep("lambda2_k",names(vast.output$Opt$par))])
	}

	if(length(grep(":",colnames(Q_ik)))>0)
	{
		# if there is an interaction, calculate factor level effects & map columns to particular factors
		# condense to a matrix with as many columns as factors
		coef.mat.1 = matrix(0,nrow=nrow(Q_ik),ncol=length(coef.names))
		coef.mat.2 = matrix(0,nrow=nrow(Q_ik),ncol=length(coef.names))
		colnames(coef.mat.2) = colnames(coef.mat.1) = coef.names

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
				coef.mat.1[,coef.names[i]] = coef.mat.1[,coef.names[i]] + rowSums(Q_ik[,cols.tmp]) * par1.tmp[j]
				coef.mat.2[,coef.names[i]] = coef.mat.2[,coef.names[i]] + rowSums(Q_ik[,cols.tmp]) * par2.tmp[j]
				rm(list=c("cols.tmp"))
			}
			rm(list=c("names.tmp","lvls.tmp","par1.tmp","par2.tmp"))
		}

		# append column for interaction effect
			coef.mat.1 = cbind(coef.mat.1,rowSums(t(apply(Q_ik[,grep(":",colnames(Q_ik))],1,function(x)as.vector(x)*par.1[grep(":",colnames(Q_ik))]))))
			coef.mat.2 = cbind(coef.mat.2,rowSums(t(apply(Q_ik[,grep(":",colnames(Q_ik))],1,function(x)as.vector(x)*par.2[grep(":",colnames(Q_ik))]))))
			coef.names = c(coef.names,"Interaction")
			level.names = c(level.names,"Interaction")
			colnames(coef.mat.2) = colnames(coef.mat.1) = coef.names

		# correct for removal of first column for identifiability
			coef.mat.2[which(rowSums(Q_ik)==0),] = coef.mat.1[which(rowSums(Q_ik)==0),] = rep(NA,length(coef.names))

	} else {

		Q_ik.1 = t(apply(Q_ik,1,function(x)as.vector(x)*par.1))
		colnames(Q_ik.1) =colnames(Q_ik)

		Q_ik.2 = t(apply(Q_ik,1,function(x)as.vector(x)*par.2))
		colnames(Q_ik.2) = colnames(Q_ik)

		# coef influence matrix
		coef.mat.1 = matrix(NA,nrow=nrow(Q_ik),ncol=length(coef.names))
		coef.mat.2 = matrix(NA,nrow=nrow(Q_ik),ncol=length(coef.names))
		colnames(coef.mat.2) = colnames(coef.mat.1) = coef.names
		for(j in 1:length(coef.names))
		{	
			if(length(grep(coef.names[j],colnames(Q_ik)))==1)
			{
				coef.mat.1[,j] = Q_ik.1[,grep(coef.names[j],colnames(Q_ik))]
				coef.mat.2[,j] = Q_ik.2[,grep(coef.names[j],colnames(Q_ik))]
			} else {
				coef.mat.1[,j] = rowSums(Q_ik.1[,grep(coef.names[j],colnames(Q_ik))])
				coef.mat.2[,j] = rowSums(Q_ik.2[,grep(coef.names[j],colnames(Q_ik))])	
			}
		}


	}
	
	# mean standardize the effects
		coef.mat.1 = apply(coef.mat.1,2,function(x)x-mean(x,na.rm=TRUE))
		coef.mat.2 = apply(coef.mat.2,2,function(x)x-mean(x,na.rm=TRUE))

	# create data.table to store the different results
		dg.1 = data.table::as.data.table(data.frame(ts = TmbData$t_iz+1))
		dg.1 = cbind(dg.1,coef.mat.1)
		dg.2 = data.table::as.data.table(data.frame(ts = TmbData$t_iz+1))
		dg.2 = cbind(dg.2,coef.mat.2)
		yr.qtr.seq = seq(from=model.start.year,to=2030,by=0.25)[1:max(dg.1$ts)]

		dg.1$yr.qtr = yr.qtr.seq[dg.1$ts]
		dg.2$yr.qtr = yr.qtr.seq[dg.2$ts]

		coef.col.names = coef.names


	# check to see if vessel effect is available
		if(length(unique(TmbData$v_i+1))>1)
		{
			vessel.index = TmbData$v_i+1
			ve.1 = Report$eta1_vc
			ve.2 = Report$eta2_vc

			dg.1$vessel = ve.1[vessel.index]
			dg.2$vessel = ve.2[vessel.index]

			coef.col.names = c(coef.col.names,"vessel")
			level.names = c(level.names,"vessel")
		}



		dg.1 = dg.1[,lapply(.SD,mean,na.rm=TRUE),by=yr.qtr,.SDcols=coef.col.names][order(yr.qtr)]
		dg.1$Component = "Encounter probability"
		dg.2 = dg.2[,lapply(.SD,mean,na.rm=TRUE),by=yr.qtr,.SDcols=coef.col.names][order(yr.qtr)]
		dg.2$Component = "Positive catch"

		g = rbind(dg.1,dg.2) %>% data.table::melt(.,id.vars = c("yr.qtr","Component")) %>% .[,variable:=factor(as.character(variable),levels = coef.col.names,labels=level.names)] %>%
			ggplot2::ggplot() + ggthemes::theme_few() + ggplot2::geom_hline(yintercept = 0,color="gray70") +
			ggplot2::xlab("Year") +
			ggplot2::ylab("Influence") +
			ggplot2::ggtitle("Influence on standardized index") +
			ggplot2::geom_point(ggplot2::aes(x=yr.qtr,y=value),alpha=pt.alpha) +
			ggplot2::geom_smooth(ggplot2::aes(x=yr.qtr,y=value),se=FALSE,span=smooth.span) +
			ggplot2::facet_grid(variable~Component,drop=FALSE)
	
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

