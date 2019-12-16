


make_extrapolation_info.ndd = function(Region,strata.limits,observations_LL,input.grid,crs.en,crs.ll,strata.sp)
# modified from Jim Thorson's FishStatUtils::Prepare_User_Extrapolation_Data_Fn.R
# returns a tagged list used in other functions
# a_el
	# The area associated with each extrapolation grid cell (rows) and strata (columns)
# Data_Extrap
	# A data frame describing the extrapolation grid with columns: Lon, Lat, Area_km2, Include, E_km, N_km
	# Include appears to be set to 1 for all areas? or set it to 0 for grid cells where you don't want an extrapolation??
# Area_km2_x
	# the area associated with each row of Data_Extrap, in units square-kilometers
# maybe need to add back "zone" and "flip_around_dateline" for plotting purposes??
{
	# define input_grid
	  	# x.range = floor(range(observations_LL$Lon)/input.grid.res)*input.grid.res
	  	# y.range = floor(range(observations_LL$Lat)/input.grid.res)*input.grid.res

	  	# # x.range[1] = x.range[1] - 0.5*input.grid.res
	  	# # x.range[2] = x.range[2] + 0.5*input.grid.res
	  	# # y.range[1] = y.range[1] - 0.5*input.grid.res
	  	# # y.range[2] = y.range[2] + 0.5*input.grid.res
	  	# x.seq = seq(from=x.range[1],to=x.range[2],by=input.grid.res)
	  	# y.seq = seq(from=y.range[1],to=y.range[2],by=input.grid.res)

	  	# input_grid = expand.grid(Lon=x.seq,Lat=y.seq)
	  	# input_grid$Area_km2 = (input.grid.res*110)^2
	  	# input_grid = as.matrix(input_grid)

	# Read extrapolation data
  		Data_Extrap = input.grid

  	# Survey areas
  		Area_km2_x = Data_Extrap[,'Area_km2']
  
  	# Augment with strata for each extrapolation cell
  		if(missing(strata.sp))
  		{
  			a_el = data.frame(All_areas = Area_km2_x)
  		} else {
  			# load packages
				library(rgeos,quietly=TRUE)
			# define union of model region
				n.substrata = length(strata.sp)
				full.reg = strata.sp[1]
				for(i in 2:length(strata.sp)){full.reg = gUnion(full.reg,strata.sp[i])}
				strata.sp = rbind(full.reg,strata.sp)
				new_IDs = c("all.strata",paste0("sub.strata.",1:n.substrata))
				for (i in 1:length(slot(strata.sp, "polygons")))
				{
				  slot(slot(strata.sp, "polygons")[[i]], "ID") = new_IDs[i]
				}

				a_el = as.data.frame(matrix(NA,nrow=length(Area_km2_x),ncol=length(strata.sp)+1))
				colnames(a_el) = c("all_areas",names(strata.sp))

				extrap.points = as.data.frame(Data_Extrap)
				coordinates(extrap.points) = c("Lon", "Lat")
	  			proj4string(extrap.points) = proj4string(strata.sp)
				a_el[,1] = Area_km2_x
				for(i in 1:length(strata.sp))
				{
					a_el[,i+1] = Area_km2_x * ifelse(is.na(over(extrap.points,strata.sp[i])),0,1)
				}
  		}


  # Convert extrapolation-data to an Eastings-Northings coordinate system
  		tmpEN = Convert_LL_to_EastNorth_Fn.ndd( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'],crs.en=crs.en,crs.ll=crs.ll)                                                         #$

  # Extra junk
	  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
	  if( all(c("E_km","N_km") %in% colnames(Data_Extrap)) ){
	    Data_Extrap[,c('E_km','N_km')] = tmpEN[,c('E_km','N_km')]
	  }else{
	    Data_Extrap = cbind( Data_Extrap, 'E_km'=tmpEN[,'E_km'], 'N_km'=tmpEN[,'N_km'] )
	  }

  # Return
	  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=NA, "flip_around_dateline"=NA, "Area_km2_x"=Area_km2_x)
	  return( Return )
}