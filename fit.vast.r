

# Nicholas Ducharme-Barth
# 10/12/2019
# fit VAST model
# compatible with VAST v8_3_0 through 16/12/2019

# RunDir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/2019_SKJ_manuscript/VAST/model_runs/"
# SaveDir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/2019_SKJ_manuscript/VAST/model_runs/"
# SourceDir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/2019_SKJ_manuscript/VAST/" 
# Q_ik = NULL
# vf.re = FALSE
# enviro = 0
# FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1)
# RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)
# ObsModel_ez = c(1,3)
# fine_scale=TRUE
# input.grid.res=1
# knot.config = 1
# n_x=100
# grid_size_km=100
# Version="VAST_v8_3_0"
# Method="Mesh"
# strata.sp = skj.alt2019.shp
# enviro=enviro.a

fit.vast = function(Data_Geostat,RunDir,SaveDir,SourceDir,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=TRUE,input.grid.res=1,knot.config = 1,n_x=100,grid_size_km=100,Version="VAST_v8_3_0",Method="Mesh",strata.sp,enviro)
{
	A = proc.time()

	library(TMB)
	library(VAST)
	library(FishStatsUtils)
	library(data.table)
	library(sp)
	source(paste0(SourceDir,"ndd.VAST.utils.r"))

	if(!dir.exists(RunDir)){dir.create(RunDir, recursive = TRUE)}
		setwd(RunDir)
	if(!dir.exists(SaveDir)){dir.create(SaveDir, recursive = TRUE)}


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
			Extrapolation_List = make_extrapolation_info.ndd(Region = Region,strata.limits = strata.limits,observations_LL = Data_Geostat[,c("Lat","Lon")],input.grid = input.grid, crs.en = crs.en,crs.ll = crs.ll,strata.sp)

		# punch-out all extrapolation grid cells that are on land and outside of "data region"
			smooth.hull = smooth.hull.sp(Data_Geostat[,c("Lon","Lat")],crs.ll=crs.ll,buffer.ll=2.5,vertices = 45, k=5)
			extrap.df = as.data.frame(Extrapolation_List$Data_Extrap)
			coordinates(extrap.df) = ~E_km + N_km
			proj4string(extrap.df) = crs.en
			# plot(smooth.hull)

		# find overlap with land
		# load land shape file
			load(paste0(SourceDir,"pacific.coast.RData"))
			pacific.coast = spTransform(pacific.coast,crs.en)
			over.index.land = which(is.na(over(extrap.df,pacific.coast)))
		# find overlap with data hull
			smooth.hull.trans = spTransform(smooth.hull,crs.en)
		# identify number of extrapolation cells (within data hull and not on land) corresponding to each knot
			over.index.region = which(!is.na(over(extrap.df,smooth.hull.trans)))
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
												  iter.max=1000, randomseed=seed, nstart=100, DirPath=SaveDir, Save_Results=FALSE,crs.en = crs.en,crs.ll = crs.ll)

			Data_Geostat = cbind(Data_Geostat, knot_i = Spatial_List$knot_i)

		# plot spatial domain of the model
			MapDetails_List = make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, 
			"spatial_list"=Spatial_List, "Extrapolation_List"=Extrapolation_List )

			plot_data(Extrapolation_List, Spatial_List, Data_Geostat = Data_Geostat, PlotDir = SaveDir)

		# allow for environmental covariates
			if(missing(enviro))
	  		{
	  			X_gtp = NULL
	  			X_itp = NULL
	  		} else {
	  			library(splines)
	  			enviro.format = make_covariates.ndd(formula = as.formula(enviro[["formula"]]),covariate_data=enviro[["covariate_data"]], Year_i=Data_Geostat[,'Year'], spatial_list=Spatial_List, extrapolation_list=Extrapolation_List)
	  			X_gtp = enviro.format$X_gtp
	  			X_itp = enviro.format$X_itp
	  		}

		# run the model
			TmbData = make_data("Version"=Version, "FieldConfig"=FieldConfig, "OverdispersionConfig"=OverdispersionConfig, 
							"RhoConfig"=RhoConfig, "ObsModel_ez"=ObsModel_ez, 
							"c_iz"=as.numeric(Data_Geostat[,'Spp'])-1, "b_i"=Data_Geostat[,'Response_variable'], 
							"a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=v_i, "Q_ik" = Q_ik, "X_gtp" = X_gtp, "X_itp" = X_itp,
							"t_iz"=Data_Geostat[,'Year'], 
							"Options"=Options, "spatial_list" = Spatial_List )

			TmbList = make_model("TmbData"=TmbData, "RunDir"=RunDir, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Spatial_List$Method)

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
			save(vast_output,file=paste0(SaveDir,"vast_output.RData"))
		return(vast_output)
}
