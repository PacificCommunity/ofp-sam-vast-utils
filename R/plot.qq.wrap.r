

#' Plot Q-Q plot for aggregate model fit. Currently only available for delta-lognormal model.
#' 
#' @param vast_output Output from a call to FishStatsUtils::fit_model
#' @param error.structure Character string denoting the error structure used to fit the model. This is needed to calculate the correct quantiles for the residuals.
#' @param nsims The number of random draws to use when calculating the quantiles. Default is 1000.
#' @param save.dir Path to the directory where the outputs will be saved
#' @export
#' @importFrom scales alpha
#' @importFrom boot inv.logit



plot.qq.wrap = function(vast_output,error.structure="dln",nsims=1000,save.dir)
{	
	if(!(error.structure %in% c("dln","dg")))
	{
		stop("Only defined for models not using a delta-lognormal or delta-gamma error structure")
	}

	######################################################
	# define internal functions
		calc.qq.binom = function(vast_output)
		# modified from Laura Tremblay-Boyer code used in 2018 alb assessment 
		{

		    p = boot::inv.logit(vast_output$Report$P1_iz)
		    y = ifelse(vast_output$data_frame$b_i>0,1,0)
		    n = rep(1, length(y)) ## assuming no prior weights

		    y = n * y
		    a = pbinom(y - 1, n, p)
		    b = pbinom(y, n, p)
		    u = runif(n = length(y), min = a, max = b)
		    bino.res = qnorm(u)
		    return(bino.res)
		}

		calc.qq.pos = function (vast_output, nsim=1000,error.structure)
		# modified from Laura Tremblay-Boyer code used in 2018 alb assessment
		# Jim Thorson's FishStatsUtils::plot_quantile_diagnostic
		{
		    TmbData = vast_output$data_frame
		    Report = vast_output$Report
		    Which = which(TmbData$b_i > 0) ## index of positive sets

		    pow = function(a, b){ a^b}
		    ## set-up
		    Q = rep(NA, length(Which))
		    y = array(NA, dim = c(length(Which), nsim))
		    var_y = rep(NA, length(Which))
		    ## get prediction of indiviuals caught (rate * area swept)
		   		pred_y = exp(Report$P2_iz[Which]) * TmbData$a_i[Which]
		   		if(error.structure == "dln")
		   		{
			        mean_y = log(pred_y) - pow(Report$SigmaM[1,1],2)/2 # corrected mean is log_exp - 0.5sigma^2
			        y = rlnorm(n=nsim*length(Which), meanlog=rep(mean_y, each=nsim), sdlog=Report$SigmaM[1,1])
			        dim(y) = c(nsim, length(Which)) ## draws in rows, observations in columns
			        var_y = apply(y, 2, var) ## estimate of variance
			        stdres = (log(TmbData$b_i[Which]) - mean_y)
		   		} else if (error.structure == "dg"){
		   			# stop("I lied, I haven't implemented this yet. Try again with error.structure = 'dln'.")
		   			gamma.scale = pow(Report$SigmaM[1,1],2) * pred_y;
		   			mean_y = gamma.scale*(1/pow(sigmaM[i_e,1],2))
         			y = rgamma(n=nsim*length(Which), shape=rep(1/pow(sigmaM[i_e,1],2), each=nsim), scale=gamma.scale)
         			dim(y) = c(nsim, length(Which)) ## draws in rows, observations in columns
			        var_y = apply(y, 2, var) ## estimate of variance
			        stdres = TmbData$b_i[Which] - mean_y
		   		} else {
		   			stop("How did you sneak an other error structure into this function?? Try again with error.structure = 'dln' or 'dg'.")
		   		}

		        pos.res = list(var_y = var_y, pred_y = pred_y, stdres=stdres,nsim=nsim)
		        pos.res$which.pos = Which

		    return(pos.res)    
		}
	######################################################
	######################################################

	# plotting
		if(missing(save.dir))
		{
    			par(mfrow=c(1,2), mar=c(4,5,1,1), las=1)
	   			qq.bin = calc.qq.binom(vast_output)
	    		ab = qqnorm(qq.bin, plot.it=FALSE)
		    	qqnorm(qq.bin, main='Encounter-rate component',col=scales::alpha('dodgerblue4',0.3),cex.axis=1.5,cex.lab=1.5)
		   		q2get = c(0.001,0.01,0.025,0.975,0.99,0.999)
	            q2p = quantile(ab$x, q2get)
	            abline(v=q2p, col='azure3')
	            pu = par('usr')
	            text(q2p, pu[4]-0.15, q2get, offset=0.3, col='azure3', srt=90, pos=2, xpd=NA, cex=0.8)
		    	qqline(qq.bin, col='grey50')
		    	rm(list=c("ab","q2get","q2p","pu"))
				pos.res = calc.qq.pos(vast_output,nsims,error.structure)
				ab = qqnorm(pos.res$stdres, plot.it=FALSE)
				plot(ab$x, ab$y, type="n", col=scales::alpha('dodgerblue4',0.3), xlab='Theoretical quantiles', ylab='Sample quantiles', main='Positive catch-rate component',las=1,cex.axis=1.5,cex.lab=1.5)
	   			qqline(pos.res$stdres, col='grey50')
	            q2get = c(0.001,0.01,0.025,0.975,0.99,0.999)
	            q2p = quantile(ab$x, q2get)
	            abline(v=q2p, col='azure3')
	            pu = par('usr')
	            text(q2p, pu[4]-0.15, q2get, offset=0.3, col='azure3', srt=90, pos=2, xpd=NA, cex=0.8)
				points(ab$x, ab$y, col=scales::alpha('dodgerblue4',0.3))
		} else {
			if (! dir.exists(save.dir))dir.create(save.dir,recursive=TRUE)
			    png(filename=paste0(save.dir,"qq.",error.structure,".png"),width = 9.5, height = 5.5, res=300 ,units = "in",  bg = "white", type = "windows")
    			par(mfrow=c(1,2), mar=c(4,5,1,1), las=1)
	   			qq.bin = calc.qq.binom(vast_output)
	    		ab = qqnorm(qq.bin, plot.it=FALSE)
		    	qqnorm(qq.bin, main='Encounter-rate component',col=scales::alpha('dodgerblue4',0.3),cex.axis=1.5,cex.lab=1.5)
		   		q2get = c(0.001,0.01,0.025,0.975,0.99,0.999)
	            q2p = quantile(ab$x, q2get)
	            abline(v=q2p, col='azure3')
	            pu = par('usr')
	            text(q2p, pu[4]-0.15, q2get, offset=0.3, col='azure3', srt=90, pos=2, xpd=NA, cex=0.8)
		    	qqline(qq.bin, col='grey50')
		    	rm(list=c("ab","q2get","q2p","pu"))
				pos.res = calc.qq.pos(vast_output,nsims,error.structure)
				ab = qqnorm(pos.res$stdres, plot.it=FALSE)
				plot(ab$x, ab$y, type="n", col=scales::alpha('dodgerblue4',0.3), xlab='Theoretical quantiles', ylab='Sample quantiles', main='Positive catch-rate component',las=1,cex.axis=1.5,cex.lab=1.5)
	   			qqline(pos.res$stdres, col='grey50')
	            q2get = c(0.001,0.01,0.025,0.975,0.99,0.999)
	            q2p = quantile(ab$x, q2get)
	            abline(v=q2p, col='azure3')
	            pu = par('usr')
	            text(q2p, pu[4]-0.15, q2get, offset=0.3, col='azure3', srt=90, pos=2, xpd=NA, cex=0.8)
				points(ab$x, ab$y, col=scales::alpha('dodgerblue4',0.3))
				dev.off()
		}
}

