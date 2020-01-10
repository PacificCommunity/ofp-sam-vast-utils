

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
#' @param input.grid.res Resolution of extrapolation grid in kilometers. d
#' @param knot_method knot_method whether to determine location of GMRF vertices based on the location of samples \code{knot_method=`samples`} or extrapolation-grid cells within the specified strata \code{knot_method='grid'}
#' @param n_x Number of knots
#' @param Version Version of VAST to use. Compatible with version "VAST_v8_3_0"
#' @param Method Method to use for defining spatial field. default setting = "Mesh"
#' @param strata.sp [Optional] If present, a shapefile containing the strata boundaries to calculate the indicies for
#' @param enviro [Optional] If present, a named-list of length two is required: "formula" is a character string that can be coerced to a formula using \code{stats::as.formula}, and "covariate_data" is a data frame with the following columns - Year, Lon, Lat, covariates...
#' @return Named list "vast_output"
#' \describe{
#'	\item{Opt}{the diagnostics from the model run}
#'  \item{Report}{the objects estimated and calculated by VAST}
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
#' @importFrom FishStatsUtils plot_data
#' @importFrom VAST make_data
#' @importFrom VAST make_model
#' @importFrom TMBhelper fit_tmb


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

fit.vast = function(Data_Geostat,RunDir,SaveDir,save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=TRUE,input.grid.res=1,knot_method = "grid",n_x=100,Version="VAST_v8_3_0",Method="Mesh",strata.sp,enviro)
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
	
	# Jim Thorson uses 110 as the approximation for the number of km in a degree of lat/lon
			input.grid = as.matrix(expand.grid(Lat = seq(from=grid_bounds[1],to=grid_bounds[2],by=input.grid.res),
					Lon = seq(from=grid_bounds[3],to=grid_bounds[4],by=input.grid.res), Area_km2 = grid_size_km))
			crs.en = paste0("+proj=tpeqd +lat_1=",mean(grid_bounds[1:2])," +lon_1=",round(grid_bounds[3] + (1/3)*abs(diff(grid_bounds[3:4])))," +lat_2=",mean(grid_bounds[1:2])," +lon_2=",round(grid_bounds[3] + (2/3)*abs(diff(grid_bounds[3:4])))," +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
			crs.ll = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
			if(missing(strata.sp))
			{
				Extrapolation_List = make_extrapolation_info.ndd(Region = Region,strata.limits = strata.limits,input.grid = input.grid, crs.en = crs.en,crs.ll = crs.ll)
			} else {
				Extrapolation_List = make_extrapolation_info.ndd(Region = Region,strata.limits = strata.limits,input.grid = input.grid, crs.en = crs.en,crs.ll = crs.ll,strata.sp)
			}		 

		# punch-out all extrapolation grid cells that are on land and outside of "data region"
			smooth.hull = smooth.hull.sp(Data_Geostat[,c("Lon","Lat")],crs.ll=crs.ll,buffer.ll=2.5,d.scalar = 0.15)
			extrap.df = as.data.frame(Extrapolation_List$Data_Extrap)
			sp::coordinates(extrap.df) = ~E_km + N_km
			sp::proj4string(extrap.df) = crs.en
			# plot(smooth.hull)

		# find overlap with land
		# load land shape file
			data("pacific.coast")
			pacific.coast = sp::spTransform(pacific.coast,crs.en)
			over.index.land = which(is.na(sp::over(extrap.df,pacific.coast)))
		# find overlap with data hull
			smooth.hull.trans = sp::spTransform(smooth.hull,crs.en)
		# identify number of extrapolation cells (within data hull and not on land) corresponding to each knot
			over.index.region = which(!is.na(sp::over(extrap.df,smooth.hull.trans)))
			over.index = intersect(over.index.land,over.index.region)

		# modify Extrapolation_List
			Extrapolation_List$a_el = data.frame(All_areas = Extrapolation_List$a_el[over.index,])
			Extrapolation_List$Data_Extrap = Extrapolation_List$Data_Extrap[over.index,]
			Extrapolation_List$Area_km2_x = Extrapolation_List$Area_km2_x[over.index]

		# with the modified spatial list function be sure to pass to Lon_i and Lat_i the lat and lon transformed to N_km and E_km using Convert_LL_to_EastNorth_Fn.ndd()
			seed = 123 
			# ll_to_EN = Convert_LL_to_EastNorth_Fn.ndd( Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'],crs.en = crs.en,crs.ll = crs.ll)
			Spatial_List = make_spatial_info.ndd(n_x=n_x, Lon_i=Data_Geostat[,'Lon'], Lat_i=Data_Geostat[,'Lat'], Extrapolation_List = Extrapolation_List, knot_method="grid", Method="Mesh",
												  grid_size_km=grid_size_km, grid_size_LL=input.grid.res, fine_scale=fine_scale, Network_sz_LL=NULL,
												  iter.max=1000, randomseed=seed, nstart=100, DirPath=SaveDir, Save_Results=save.output,crs.en = crs.en,crs.ll = crs.ll)

			Data_Geostat = cbind(Data_Geostat, knot_i = Spatial_List$knot_i)

		# prep info for plotting spatial domain of the model
			MapDetails_List = FishStatsUtils::make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, 
			"spatial_list"=Spatial_List, "Extrapolation_List"=Extrapolation_List )

		# allow for environmental covariates
			if(missing(enviro))
	  		{
	  			X_gtp = NULL
	  			X_itp = NULL
	  		} else {
	  			enviro.format = make_covariates.ndd(formula = as.formula(enviro[["formula"]]),covariate_data=enviro[["covariate_data"]], Year_i=Data_Geostat[,'Year'], spatial_list=Spatial_List, extrapolation_list=Extrapolation_List)
	  			X_gtp = enviro.format$X_gtp
	  			X_itp = enviro.format$X_itp
	  		}

		# run the model
			TmbData = VAST::make_data("Version"=Version, "FieldConfig"=FieldConfig, "OverdispersionConfig"=OverdispersionConfig, 
							"RhoConfig"=RhoConfig, "ObsModel_ez"=ObsModel_ez, 
							"c_iz"=as.numeric(Data_Geostat[,'Spp'])-1, "b_i"=Data_Geostat[,'Response_variable'], 
							"a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=v_i, "Q_ik" = Q_ik, "X_gtp" = X_gtp, "X_itp" = X_itp,
							"t_iz"=Data_Geostat[,'Year'], 
							"Options"=Options, "spatial_list" = Spatial_List )

			TmbList = VAST::make_model("TmbData"=TmbData, "RunDir"=RunDir, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Spatial_List$Method)

		# Check parameters
			Obj = TmbList[["Obj"]]
			Obj$fn( Obj$par )
			Obj$gr( Obj$par )

		# Estimate fixed effects and predict random effects
			Opt = TMBhelper::fit_tmb( obj = Obj, lower = TmbList[["Lower"]], upper = TmbList[["Upper"]],
				getsd = FALSE, savedir = NULL, bias.correct = FALSE, newtonsteps = 3 )
			B = proc.time()
			Report = Obj$report()
			fit.time = (B-A)[3]

			vast_output = list( "Opt"=Opt, "Report"=Report, "TmbData"=TmbData, "Extrapolation_List"=Extrapolation_List, "fit.time"=fit.time,"MapDetails_List"=MapDetails_List )
			if(save.output)
			{
				save(vast_output,file=paste0(SaveDir,"vast_output.RData"))
			}
			setwd(origwd)

		return(vast_output)
}
