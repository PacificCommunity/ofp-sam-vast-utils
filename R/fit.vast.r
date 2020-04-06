

#' Fits a VAST model, compatible with VAST v8_3_0 through 16/12/2019
#' Some of the descriptions come from JT VAST package descriptions
#' 
#' @param Data_Geostat A data-frame of i rows containing the following columns: Response_variable, Year, Lon, Lat, Spp, AreaSwept_km2, Vessel
#' @param RunDir Path to the directory where the .cpp VAST source code is stored or compiled
#' @param SaveDir Path to the directory where the outputs will be saved
#' @param save.output TRUE or FALSE
#' @param Q_ik Matrix or i rows and k covariates impacting catchability. Can be created using \code{stats::model.matrix}
#' @param vf.re TRUE or FALSE switch indicating if vessel random effects are to be estimated. If so then the Vessel column in Data_Geostat is used
#' @param FieldConfig Controls the number of factors estimated with the spatial and spatiotemporal random fields. default setting = c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1)
#' @param RhoConfig Controls the temporal structure of the annual intercepts and the spatiotemporal random field. default setting = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)
#' @param ObsModel_ez Controls the error structure. default setting = c(1,3)
#' \describe{
#'   \item{ObsModel_ez[e,1]=0}{Normal}
#'   \item{ObsModel_ez[e,1]=1}{Lognormal}
#'   \item{ObsModel_ez[e,1]=2}{Gamma}
#'   \item{ObsModel_ez[e,1]=3}{Inverse-Gaussian}
#'   \item{ObsModel_ez[e,1]=5}{Negative binomial}
#'   \item{ObsModel_ez[e,1]=6}{Conway-Maxwell-Poisson (likely to be very slow)}
#'   \item{ObsModel_ez[e,1]=7}{Poisson (more numerically stable than negative-binomial)}
#'   \item{ObsModel_ez[e,1]=8}{Compound-Poisson-Gamma, where the expected number of individuals is the 1st-component, the expected biomass per individual is the 2nd-component, and SigmaM is the variance in positive catches (likely to be very slow)}
#'   \item{ObsModel_ez[e,1]=9}{Binned-Poisson (for use with REEF data, where 0=0 individual; 1=1 individual; 2=2:10 individuals; 3=>10 individuals)}
#'   \item{ObsModel_ez[e,1]=10}{Tweedie distribution, where epected biomass (lambda) is the product of 1st-component and 2nd-component, variance scalar (phi) is the 1st component, and logis-SigmaM is the power}
#'   \item{ObsModel_ez[e,1]=11}{Zero-inflated Poisson with additional normally-distributed variation overdispersion in the log-intensity of the Poisson distribution}
#'   \item{ObsModel_ez[e,1]=12}{Poisson distribution (not zero-inflated) with log-intensity from the 1st linear predictor, to be used in combination with the Poisson-link delta model for combining multiple data types}
#'   \item{ObsModel_ez[e,1]=13}{Bernoilli distribution using complementary log-log (cloglog) link from the 1st linear predictor, to be used in combination with the Poisson-link delta model for combining multiple data types}
#'   \item{ObsModel_ez[e,1]=14}{Similar to 12, but also including lognormal overdispersion}
#'   \item{ObsModel_ez[e,2]=0}{Conventional delta-model using logit-link for encounter probability and log-link for positive catch rates}
#'   \item{ObsModel_ez[e,2]=1}{Alternative "Poisson-link delta-model" using log-link for numbers-density and log-link for biomass per number}
#'   \item{ObsModel_ez[e,2]=2}{Link function for Tweedie distribution, necessary for \code{ObsModel_ez[e,1]=8} or \code{ObsModel_ez[e,1]=10}}
#'   \item{ObsModel_ez[e,2]=3}{Conventional delta-model, but fixing encounter probability=1 for any year where all samples encounter the species}
#'   \item{ObsModel_ez[e,2]=4}{Poisson-link delta-model, but fixing encounter probability=1 for any year where all samples encounter the species and encounter probability=0 for any year where no samples encounter the species}
#' }
#' @param fine_scale TRUE or FALSE. Better maps and slightly better index fit when TRUE but is slower.
#' @param input.grid.res Resolution of extrapolation grid in kilometers.
#' @param crop.extrap.by.data TRUE or FALSE: If TRUE then the extrapolation grid is cropped by the smooth hull surrounding the data
#' @param knot_method knot_method whether to determine location of GMRF vertices based on the location of samples \code{knot_method=`samples`} or extrapolation-grid cells within the specified strata \code{knot_method='grid'}
#' @param n_x Number of knots
#' @param Version Version of VAST to use. Compatible with version "VAST_v8_3_0"
#' @param Method Method to use for defining spatial field. default setting = "Mesh"
#' @param ADREPORT TRUE or FALSE. Calculate the SD for the params and index?
#' @param normalize_idx TRUE or FALSE. Normalize the index (and the SE) by dividing by the mean of the index
#' @param Xconfig_zcp OPTIONAL, 3D array of settings for each dynamic density covariate, where the first dimension corresponds to 1st or 2nd linear predictors, second dimension corresponds to model category, and third dimension corresponds to each density covariate
#' \describe{
#'   \item{Xconfig_zcp[z,c,p]=0}{\code{X_itp[,,p]} has no effect on linear predictor z for category c}
#'   \item{Xconfig_zcp[z,c,p]=1}{\code{X_itp[,,p]} has a linear effect on linear predictor z for category c}
#'   \item{Xconfig_zcp[z,c,p]=2}{\code{X_itp[,,p]} has a spatially varying, zero-centered linear effect on linear predictor z for category c}
#'   \item{Xconfig_zcp[z,c,p]=3}{\code{X_itp[,,p]} has a spatially varying linear effect on linear predictor z for category c}
#' }
#' @param slim.output TRUE or FALSE, if true then vast_output only contains idx and/or idx.se, fit.time, mgc
#' @param strata.sp [Optional] If present, a shapefile containing the strata boundaries to calculate the indicies for
#' @param enviro [Optional] If present, a named-list of length two is required: "formula" is a character string that can be coerced to a formula using \code{stats::as.formula}, and "covariate_data" is a data frame with the following columns - Year, Lon, Lat, covariates...
#' @return Named list "vast_output"
#' \describe{
#'  \item{idx}{the index within each strata if strata is provided}
#'  \item{idx.se}{the associated se for the index}
#'	\item{Opt}{the diagnostics from the model run}
#'  \item{Report}{the objects estimated and calculated by VAST}
#' 	\item{Sdreport}{the report generated by ADREPORT}
#'  \item{TmbData}{the data passed to TMB and used to fit the model}
#'  \item{Extrapolation_List}{the extrapolation list}
#'  \item{fit.time}{the time to run the function}
#'  \item{MapDetails_List}{the map details to make additional plots}
#' }
#' @export
#' @import splines
#' @importFrom sp coordinates
#' @importFrom sp proj4string
#' @importFrom sp over
#' @importFrom sp spTransform
#' @importFrom FishStatsUtils make_map_info
#' @importFrom FishStatsUtils make_covariates
#' @importFrom FishStatsUtils plot_data
#' @importFrom FishStatsUtils make_extrapolation_info
#' @importFrom FishStatsUtils make_spatial_info
#' @importFrom VAST make_data
#' @importFrom VAST make_model
#' @importFrom TMBhelper fit_tmb
#' @importFrom TMB sdreport
#' @importFrom TMB summary.sdreport


### Function defaults for testing 
# RunDir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/2019_SKJ_manuscript/VAST/model_runs/"
# SaveDir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/2019_SKJ_manuscript/VAST/model_runs/"
# SourceDir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/2019_SKJ_manuscript/VAST/" 
# Q_ik = NULL
# vf.re = FALSE
# FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1)
# RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)
# ObsModel_ez = c(1,3)
# fine_scale=TRUE
# input.grid.res=1
# knot_method = "grid"
# n_x=100
# Version="VAST_v8_3_0"
# Method="Mesh"
# strata.sp = skj.alt2019.shp
# enviro=enviro.a

fit.vast = function(Data_Geostat,RunDir,SaveDir,save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=TRUE,input.grid.res=1,crop.extrap.by.data=TRUE,knot_method = "grid",n_x=100,Version="VAST_v8_3_0",Method="Mesh",ADREPORT=TRUE,normalize_idx=FALSE,Xconfig_zcp=NULL,slim.output=FALSE,strata.sp,enviro, newton_steps = 3)
{
	A = proc.time()

	if(!dir.exists(RunDir)){dir.create(RunDir, recursive = TRUE)}
		origwd = getwd()
		setwd(RunDir)
	if(save.output)
	{
		if(!dir.exists(SaveDir)){dir.create(SaveDir, recursive = TRUE)}
	}


	# Decide which post-hoc calculations to include in VAST output
		Options = c('SD_site_density'=FALSE, 'SD_site_logdensity'=FALSE, 'Calculate_Range'=FALSE, 'SD_observation_density'=FALSE, 'Calculate_effective_area'=FALSE,
    'Calculate_Cov_SE'=FALSE, 'Calculate_Synchrony'=FALSE, 'Calculate_Coherence'=FALSE, 'Calculate_proportion'=FALSE, 'normalize_GMRF_in_CPP'=TRUE,
    'Calculate_Fratio'=FALSE, 'Estimate_B0'=FALSE, 'Project_factors'=FALSE, 'treat_nonencounter_as_zero'=FALSE, 'simulate_random_effects'=TRUE )
		Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )
		if(vf.re)
		{
			OverdispersionConfig = c("Eta1"=1, "Eta2"=1)
			v_i = Data_Geostat$Vessel
		} else {
			OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
			v_i = NULL
		}

	# Determine strata within the study region 
			strata.limits = data.frame('STRATA'="All_areas")
			Region = "user"

	# Set spatial extent and extrapolation grid
			grid_size_km = (110 * input.grid.res)^2 # the distance between grid cells for the 2D AR1 grid if Method == "Grid"
			grid_bounds = c(floor(min(Data_Geostat[,c("Lat")])/5)*5,ceiling(max(Data_Geostat[,c("Lat")])/5)*5,floor(min(Data_Geostat[,c("Lon")])/5)*5,ceiling(max(Data_Geostat[,c("Lon")])/5)*5)
			crs.en = paste0("+proj=tpeqd +lat_1=",mean(grid_bounds[1:2])," +lon_1=",round(grid_bounds[3] + (1/3)*abs(diff(grid_bounds[3:4])))," +lat_2=",mean(grid_bounds[1:2])," +lon_2=",round(grid_bounds[3] + (2/3)*abs(diff(grid_bounds[3:4])))," +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
			crs.ll = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

	# Jim Thorson uses 110 as the approximation for the number of km in a degree of lat/lon
			tmp.input.grid = expand.grid(Lat = seq(from=grid_bounds[1],to=grid_bounds[2],by=input.grid.res),
					Lon = seq(from=grid_bounds[3],to=grid_bounds[4],by=input.grid.res))
			input.grid = as.matrix(data.frame(Lat = tmp.input.grid$Lat,Lon = tmp.input.grid$Lon, Area_km2=area.cell(lond=tmp.input.grid$Lon,latd=tmp.input.grid$Lat,cell.size=input.grid.res,crs.ll=crs.ll)))	
			
			Extrapolation_List = FishStatsUtils::make_extrapolation_info( Region, projargs=crs.en, strata.limits=data.frame('STRATA'="All_areas"),
  															create_strata_per_region=FALSE, input_grid=input.grid, observations_LL=Data_Geostat[,c("Lat","Lon")], 
  															grid_dim_km=c(grid_size_km,grid_size_km),flip_around_dateline = TRUE)

				# reformat Extrapolation_List$a_el
					if(!missing(strata.sp))
					{
						# define union of model region
						n.substrata = length(strata.sp)
						full.reg = strata.sp[1]
						for(i in 2:length(strata.sp)){full.reg = rgeos::gUnion(full.reg,strata.sp[i])}
						strata.sp = rbind(full.reg,strata.sp)
						new_IDs = c("all.strata",paste0("sub.strata.",1:n.substrata))
						for (i in 1:length(slot(strata.sp, "polygons")))
						{
						  slot(slot(strata.sp, "polygons")[[i]], "ID") = new_IDs[i]
						}

						Extrapolation_List$a_el = as.data.frame(matrix(NA,nrow=length(Extrapolation_List$Area_km2_x),ncol=length(strata.sp)+1))
						colnames(Extrapolation_List$a_el) = c("all_areas",names(strata.sp))

						extrap.points = as.data.frame(input.grid)
						sp::coordinates(extrap.points) = c("Lon", "Lat")
			  			sp::proj4string(extrap.points) = sp::proj4string(strata.sp)
						Extrapolation_List$a_el[,1] = Extrapolation_List$Area_km2_x
						for(i in 1:length(strata.sp))
						{
							Extrapolation_List$a_el[,i+1] = Extrapolation_List$Area_km2_x * ifelse(is.na(sp::over(extrap.points,strata.sp[i])),0,1)
						} 
					}
								 

		# punch-out all extrapolation grid cells that are on land and outside of "data region"
			extrap.df = as.data.frame(Extrapolation_List$Data_Extrap)
			sp::coordinates(extrap.df) = ~E_km + N_km
			sp::proj4string(extrap.df) = crs.en

		# find overlap with land
		# load land shape file
			data("pacific.coast")
			pacific.coast = sp::spTransform(pacific.coast,crs.en)
			over.index.land = which(is.na(sp::over(extrap.df,pacific.coast)))
		if(crop.extrap.by.data)
		{
			# find overlap with data hull
			smooth.hull = smooth.hull.sp(Data_Geostat[,c("Lon","Lat")],crs.ll=crs.ll,buffer.ll=2.5,d.scalar = 0.15)
			# plot(smooth.hull)
			smooth.hull.trans = sp::spTransform(smooth.hull,crs.en)
			# identify number of extrapolation cells (within data hull and not on land) corresponding to each knot
				over.index.region = which(!is.na(sp::over(extrap.df,smooth.hull.trans)))
				over.index = intersect(over.index.land,over.index.region)
		} else {
			over.index = over.index.land
		}
		
		# modify Extrapolation_List
			Extrapolation_List$a_el = data.frame(All_areas = Extrapolation_List$a_el[over.index,])
			Extrapolation_List$Data_Extrap = Extrapolation_List$Data_Extrap[over.index,]
			Extrapolation_List$Area_km2_x = Extrapolation_List$Area_km2_x[over.index]

		# with the modified spatial list function be sure to pass to Lon_i and Lat_i the lat and lon transformed to N_km and E_km using Convert_LL_to_EastNorth_Fn.ndd()
			seed = 123 
			# ll_to_EN = Convert_LL_to_EastNorth_Fn.ndd( Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'],crs.en = crs.en,crs.ll = crs.ll)
			Spatial_List = FishStatsUtils::make_spatial_info(n_x=n_x, Lon_i=Data_Geostat[,'Lon'], Lat_i=Data_Geostat[,'Lat'], Extrapolation_List = Extrapolation_List, knot_method="grid", Method="Mesh",
												  grid_size_km=grid_size_km, grid_size_LL=input.grid.res, fine_scale=fine_scale, Network_sz_LL=NULL,
												  iter.max=1000, randomseed=seed, nstart=100, DirPath=SaveDir, Save_Results=save.output)

			Data_Geostat = cbind(Data_Geostat, knot_i = Spatial_List$knot_i)

		# prep info for plotting spatial domain of the model
			MapDetails_List = FishStatsUtils::make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, 
			"spatial_list"=Spatial_List, "Extrapolation_List"=Extrapolation_List )

		# allow for environmental covariates
			if(missing(enviro))
	  		{
	  			X_gtp = NULL
	  			X_itp = NULL
	  			Xconfig_zcp = NULL
	  		} else {
	  			enviro.format = FishStatsUtils::make_covariates(formula = as.formula(enviro[["formula"]]),covariate_data=enviro[["covariate_data"]], Year_i=Data_Geostat[,'Year'], spatial_list=Spatial_List, extrapolation_list=Extrapolation_List)
	  			X_gtp = enviro.format$X_gtp
	  			X_itp = enviro.format$X_itp
	  		}

		# run the model
			TmbData = VAST::make_data("Version"=Version, "FieldConfig"=FieldConfig, "OverdispersionConfig"=OverdispersionConfig, 
							"RhoConfig"=RhoConfig, "ObsModel_ez"=ObsModel_ez, 
							"c_iz"=as.numeric(Data_Geostat[,'Spp'])-1, "b_i"=Data_Geostat[,'Response_variable'], 
							"a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=v_i, "Q_ik" = Q_ik, "X_gtp" = X_gtp, "X_itp" = X_itp, "Xconfig_zcp"=Xconfig_zcp,
							"t_iz"=Data_Geostat[,'Year'], 
							"Options"=Options, "spatial_list" = Spatial_List )

			TmbList = VAST::make_model("TmbData"=TmbData, "RunDir"=RunDir, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Spatial_List$Method)

		# Check parameters
			Obj = TmbList[["Obj"]]
			Obj$fn( Obj$par )
			Obj$gr( Obj$par )

		# Estimate fixed effects and predict random effects
			if(ADREPORT)
			{
				Opt = TMBhelper::fit_tmb( obj = Obj, lower = TmbList[["Lower"]], upper = TmbList[["Upper"]],
				getsd = TRUE, savedir = NULL, bias.correct = FALSE, newtonsteps = newton_steps )
				Report = Obj$report()
				Sdreport = TMB::sdreport(Obj)
				B = proc.time()
				fit.time = (B-A)[3]
				idx = Report$Index_cyl[1,,]
				idx.se = array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))=="Index_cyl"),2], dim=c(dim(Report$Index_cyl)), dimnames=list(NULL,NULL,NULL) )[1,,]

				if(normalize_idx)
				{
					for(j in 1:ncol(idx))
					{
						j.mean = mean(idx[,j])
						idx[,j] = idx[,j]/j.mean
						idx.se[,j] = idx.se[,j]/j.mean
					}
				}
				vast_output = list("idx"=idx,"idx.se"=idx.se, "Opt"=Opt, "Report"=Report,"Sdreport" = Sdreport, "TmbData"=TmbData,"Spatial_List"=Spatial_List, "Extrapolation_List"=Extrapolation_List, "fit.time"=fit.time,"MapDetails_List"=MapDetails_List )
		
			} else {
				Opt = TMBhelper::fit_tmb( obj = Obj, lower = TmbList[["Lower"]], upper = TmbList[["Upper"]],
				getsd = FALSE, savedir = NULL, bias.correct = FALSE, newtonsteps = newton_steps )
				Report = Obj$report()
				B = proc.time()
				fit.time = (B-A)[3]
				idx = Report$Index_cyl[1,,]
				if(normalize_idx)
				{
					idx = apply(idx,2,function(x)x/mean(x))
				}
				vast_output = list("idx"=idx, "Opt"=Opt, "Report"=Report, "TmbData"=TmbData,"Spatial_List"=Spatial_List, "Extrapolation_List"=Extrapolation_List, "fit.time"=fit.time,"MapDetails_List"=MapDetails_List )
		
			}

			if(slim.output)
			{
				vast_output = vast_output[names(vast_output) %in% c("idx","idx.se","fit.time")]
				vast_output = c(vast_output,Opt)
				vast_output$X_gtp = X_gtp
			}

			if(save.output)
			{
				save(vast_output,file=paste0(SaveDir,"vast_output.RData"))
			}
			setwd(origwd)

		return(vast_output)
}
