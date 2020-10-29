

#' Format the estimated index for use in MFCL
#' 
#' @param vast_output Output from a call to fit.vast where slim output is FALSE
#' @param agg.years Years used to define the average period when rescaling
#' @param ts.vec a vector denoting the year-quarter of each ts for the fit.vast model output
#' @param region.idx The fit.vast model can calculate the abundance trends over many regions, some of which are superfluous. Specify the column index for the regions that you care about.
#' @param region.names What should these regions be called? A character vector of names please
#' @param mean.cv When rescaling the CV for MFCL, set the mean cv to rescale to.
#' @param missing If the cv is missing for some reason (a hold over from the conventional delta-glm) set the penalty weight to this value. The smaller the value, the smaller the impact of this data point on the likelihood in MFCL.
#' @param save.dir Path to the directory where the outputs will be saved
#' @param save.name Name stem for the output, useful when saving many model outputs in the same directory
#' @export
#' @importFrom data.table as.data.table


vast.frq.index = function(vast_output,agg.years=NULL,ts.vec=seq(from=1952,to=2018.75,by=0.25),region.idx=3:11,region.names=paste0("R",1:length(region.idx)),mean.cv=0.2,missing=0.05,save.dir,save.name)
{
	if(is.null(agg.years))
	{
		# 1) calculate regional weights
			extrap.info = vast_output$Extrapolation_List$a_el[,region.idx]
			colnames(extrap.info) = region.names
			extrap.info$knot = vast_output$Spatial_List$NN_Extrap$nn.idx
			extrap.info = data.table::as.data.table(extrap.info)

			D_yx = t(vast_output$Report$D_gcy[,1,])

			# create matrix of density in each region over time (area for each knot in each region times the density at the knot and time)
			reg.wt = matrix(NA,nrow=length(ts.vec),ncol=length(region.names))
			colnames(reg.wt) = region.names

				for(j in 1:length(region.names))
				{
					tmp.dt = extrap.info[,c(region.names[j],"knot"),with=FALSE] 
					colnames(tmp.dt) = c("reg","knot")
					reg.knot.area = as.matrix(tmp.dt[,.(Area_km2_x=sum(reg)),by=knot][order(knot)])

					for(i in 1:length(ts.vec))
					{
						reg.wt[i,j] = sum(D_yx[i,] * reg.knot.area[,"Area_km2_x"])
					}

					rm(list=c("tmp.dt","reg.knot.area"))
				}

			# calculate regional weight as average density over the whole model period
			reg.wt = colMeans(reg.wt)
			reg.wt = reg.wt/sum(reg.wt)
			names(reg.wt) = region.names

		# 2) format idx and cv
			idx = vast_output$idx[,region.idx]
			se = vast_output$idx.se[,region.idx]
			cv = se/idx
			colnames(cv) = paste0(region.names,"cv")

			reg.mean = apply(idx,2,function(x)mean(x,na.rm=TRUE))
			idx = do.call("cbind",lapply(1:length(region.idx),function(i) idx[,i]/reg.mean[i]))
			colnames(idx) = region.names

			idx.std = cbind(idx,cv)
			idx.std = cbind(ts.vec,idx.std)
			colnames(idx.std)[1] = "yrqtr"
			rm(list=c("idx","se","cv"))
		# 3) format idx and penalty weights for frq
			idx = vast_output$idx[,region.idx]
			se = vast_output$idx.se[,region.idx]
			cv = se/idx
			colnames(cv) = paste0(region.names,"cv")

			global.mean = mean(idx,na.rm=TRUE)
			idx = do.call("cbind",lapply(1:length(region.idx),function(i) idx[,i]/global.mean))
			colnames(idx) = region.names

			idx.frq = cbind(idx,cv)
			idx.frq = cbind(ts.vec,idx.frq)
			colnames(idx.frq)[1] = "yrqtr"
			cv.cols = seq(from=length(region.idx)+2,by=1,length.out=length(region.idx))
			colnames(idx.frq)[cv.cols] = paste0("R",1:length(region.idx),"penwt")

			# Rescale CVs by dividing by the mean over the model period and multiplying by mean.cv
			# Convert CVs to penalty weights
			global.cv = mean(cv,na.rm=TRUE)
			idx.frq[,cv.cols] = apply(idx.std[,cv.cols],2,function(x)(x/global.cv)*mean.cv)
			for(i in 1:length(cv.cols))
			{
			    idx.frq[,cv.cols[i]] = round(1/(2 * (idx.frq[,cv.cols[i]]^2)),3)
			    idx.frq[,cv.cols[i]] = ifelse(is.infinite(idx.frq[,cv.cols[i]]),missing,idx.frq[,cv.cols[i]])
				idx.frq[,cv.cols[i]] = ifelse(idx.frq[,cv.cols[i]]<1,missing,idx.frq[,cv.cols[i]])
			}

	} else {
		# 1) calculate regional weights
			extrap.info = vast_output$Extrapolation_List$a_el[,region.idx]
			colnames(extrap.info) = region.names
			extrap.info$knot = vast_output$Spatial_List$NN_Extrap$nn.idx
			extrap.info = data.table::as.data.table(extrap.info)

			D_yx = t(vast_output$Report$D_gcy[,1,])

			# create matrix of density in each region over time (area for each knot in each region times the density at the knot and time)
			reg.wt = matrix(NA,nrow=length(ts.vec),ncol=length(region.names))
			colnames(reg.wt) = region.names

				for(j in 1:length(region.names))
				{
					tmp.dt = extrap.info[,c(region.names[j],"knot"),with=FALSE] 
					colnames(tmp.dt) = c("reg","knot")
					reg.knot.area = as.matrix(tmp.dt[,.(Area_km2_x=sum(reg)),by=knot][order(knot)])

					for(i in 1:length(ts.vec))
					{
						reg.wt[i,j] = sum(D_yx[i,] * reg.knot.area[,"Area_km2_x"])
					}

					rm(list=c("tmp.dt","reg.knot.area"))
				}

			agg.idx = which(floor(ts.vec) %in% agg.years)
			# calculate regional weight as average density over the agg.years period
			reg.wt = colMeans(reg.wt[agg.idx,])
			reg.wt = reg.wt/sum(reg.wt)
			names(reg.wt) = region.names

		# 2) format idx and cv
			idx = vast_output$idx[,region.idx]
			se = vast_output$idx.se[,region.idx]
			cv = se/idx
			colnames(cv) = paste0(region.names,"cv")

			agg.mean = apply(idx,2,function(x)mean(x[agg.idx],na.rm=TRUE))
			idx = do.call("cbind",lapply(1:length(region.idx),function(i) idx[,i]/agg.mean[i]))
			colnames(idx) = region.names

			idx.std = cbind(idx,cv)
			idx.std = cbind(ts.vec,idx.std)
			colnames(idx.std)[1] = "yrqtr"

		# 3) format idx and penalty weights for frq
			idx.frq = idx.std
			cv.cols = seq(from=length(region.idx)+2,by=1,length.out=length(region.idx))
			colnames(idx.frq)[cv.cols] = paste0("R",1:length(region.idx),"penwt")
			# Rescale indices by multiplying by regional weight
			for(j in 1:length(region.idx))
			{
				idx.frq[,j+1] = idx.frq[,j+1] * reg.wt[j] 
			}

			# Rescale CVs by dividing by the mean over the average period and multiplying by mean.cv
			# Convert CVs to penalty weights

			idx.frq[,cv.cols] = apply(idx.std[,cv.cols],2,function(x)(x/mean(x[agg.idx],na.rm=TRUE))*mean.cv)
			for(i in 1:length(cv.cols))
			{
			    idx.frq[,cv.cols[i]] = round(1/(2 * (idx.frq[,cv.cols[i]]^2)),3)
			    idx.frq[,cv.cols[i]] = ifelse(is.infinite(idx.frq[,cv.cols[i]]),missing,idx.frq[,cv.cols[i]])
				idx.frq[,cv.cols[i]] = ifelse(idx.frq[,cv.cols[i]]<1,missing,idx.frq[,cv.cols[i]])
			}

	}


	# write.out
		if(!missing(save.dir))
		{
			if(missing(save.name))
			{
				stop("How can you save the output if you haven't specified the directory? Please specify save.dir.")
			} else {
				write.csv(reg.wt,file=paste0(save.dir,save.name,"-reg.wt.csv"))
				write.csv(idx.std,file=paste0(save.dir,save.name,"-idx.std.csv"))
				write.csv(idx.frq,file=paste0(save.dir,save.name,"-idx.frq.csv"))
			}
		} 

		idx.list = list(reg.wt=reg.wt,idx.std=idx.std,idx.frq=idx.frq)
		return(idx.list)
}