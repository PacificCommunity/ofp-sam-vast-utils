
#' modified from Jim Thorson's FishStatUtils::Prepare_User_Extrapolation_Data_Fn.R
#' uses my \code{Convert_EN_to_LL_Fn.ndd} and \code{Convert_LL_to_EastNorth_Fn.ndd} functions
#' also modifies \code{strata.limits} to match \code{strata.sp} if it is present 
#' @param Region user
#' @param strata.limits data frame with strata limits
#' @param input.grid user created extrapolation grid
#' @param crs.en Character string of the crs for the E-N projection
#' @param crs.ll Character string of the crs for the current lat-lon projections
#' @param strata.sp [Optional] If present, a shapefile containing the strata boundaries to calculate the indicies for
#' @return returns a tagged list used in other functions
#' \describe{
#'   \item{a_el}{The area associated with each extrapolation grid cell (rows) and strata (columns)}
#'   \item{Data_Extrap}{A data frame describing the extrapolation grid with columns: Lon, Lat, Area_km2, Include, E_km, N_km. Include appears to be set to 1 for all areas? or set it to 0 for grid cells where you don't want an extrapolation??}
#'   \item{Area_km2_x}{the area associated with each row of Data_Extrap, in units square-kilometers}
#'   \item{zone}{}
#'   \item{flip_around_dateline}{}
#' }
#' @importFrom rgeos gUnion
#' @importFrom sp coordinates
#' @importFrom sp proj4string
#' @importFrom sp over
#' @export

make_extrapolation_info.ndd = function(Region,strata.limits,input.grid,crs.en,crs.ll,strata.sp)
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

				a_el = as.data.frame(matrix(NA,nrow=length(Area_km2_x),ncol=length(strata.sp)+1))
				colnames(a_el) = c("all_areas",names(strata.sp))

				extrap.points = as.data.frame(Data_Extrap)
				sp::coordinates(extrap.points) = c("Lon", "Lat")
	  			sp::proj4string(extrap.points) = sp::proj4string(strata.sp)
				a_el[,1] = Area_km2_x
				for(i in 1:length(strata.sp))
				{
					a_el[,i+1] = Area_km2_x * ifelse(is.na(sp::over(extrap.points,strata.sp[i])),0,1)
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