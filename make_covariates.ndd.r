

make_covariates.ndd = function( formula=~0, covariate_data, Year_i, spatial_list, extrapolation_list ){
# added `colnames(DF_ip) = c(names(sample_data),covariate_names)` to correct column mismatch
  # Errors
  if( !is.data.frame(covariate_data) ) stop("Please ensure that `covariate_data` is a data frame")
  if( !all(c("Lat","Lon","Year") %in% names(covariate_data)) ){
    stop( "`data` in `make_covariates(.)` must include columns `Lat`, `Lon`, and `Year`" )
  }

  # transform data inputs
  sample_data = data.frame( "Year"=Year_i, "Lat"=spatial_list$latlon_i[,'Lat'], "Lon"=spatial_list$latlon_i[,'Lon'] )
  covariate_names = setdiff( names(covariate_data), names(sample_data) )

  # set of years needed
  Year_Set = min(Year_i):max(Year_i)

  # extract latitude and longitude for extrapolation grid
  latlon_g = spatial_list$latlon_g

  # Create data frame of necessary size
  DF_zp = NULL
  DF_ip = cbind( sample_data, covariate_data[rep(1,nrow(sample_data)),covariate_names] )
  colnames(DF_ip) = c(names(sample_data),covariate_names)
  DF_ip[,covariate_names] = NA

  # Loop through data and extrapolation-grid
  for( tI in seq_along(Year_Set) ){

    # Subset to same year
    tmp_covariate_data = covariate_data[ which(Year_Set[tI]==covariate_data[,'Year'] | is.na(covariate_data[,'Year'])), , drop=FALSE]
    if( nrow(tmp_covariate_data)==0 ){
      stop("Year ", Year_Set[tI], " not found in `covariate_data` please specify covariate values for all years" )
    }
    #
    Which = which(Year_Set[tI]==sample_data[,'Year'])
    # Do nearest neighbors to define covariates for observations, skipping years without observations
    if( length(Which) > 0 ){
      NN = RANN::nn2( data=tmp_covariate_data[,c("Lat","Lon")], query=sample_data[Which,c("Lat","Lon")], k=1 )
      # Add to data-frame
      nearest_covariates = tmp_covariate_data[ NN$nn.idx[,1], covariate_names, drop=FALSE ]
      DF_ip[Which, covariate_names] = nearest_covariates
    }

    # Do nearest neighbors to define covariates for extrapolation grid, including years without observations
    NN = RANN::nn2( data=tmp_covariate_data[,c("Lat","Lon")], query=latlon_g[,c("Lat","Lon")], k=1 )
    # Add rows
    nearest_covariates = tmp_covariate_data[ NN$nn.idx[,1], covariate_names, drop=FALSE ]
    newrows = cbind("Year"=Year_Set[tI], latlon_g, nearest_covariates )
    DF_zp = rbind( DF_zp, newrows )
  }
  if( any(is.na(DF_ip)) ) stop("Problem with `DF_ip` in `make_covariates(.)")

  # Convert to dimensions requested
  DF = rbind( DF_ip, DF_zp )
  X = model.matrix( update.formula(formula, ~.+0), data=DF )[,,drop=FALSE]

  # Make X_ip
  X_ip = X[ 1:nrow(DF_ip), , drop=FALSE ]
  X_itp = aperm( X_ip %o% rep(1,length(Year_Set)), perm=c(1,3,2) )

  # Make X_gpt and then permute dimensions
  X_gpt = NULL
  indices = nrow(X_ip)
  for( tI in seq_along(Year_Set) ){
    indices = max(indices) + 1:nrow(latlon_g)
    if( max(indices)>nrow(X) ) stop("Check problem in `make_covariates`")
    X_gpt = abind::abind( X_gpt, X[ indices, , drop=FALSE ], along=3 )
  }
  X_gtp = aperm( X_gpt, perm=c(1,3,2) )

  # warnings
  if( any(apply(X_gtp, MARGIN=2:3, FUN=sd)>10 | apply(X_itp, MARGIN=2:3, FUN=sd)>10) ){
    warning("The package author recommends that you rescale covariates in `covariate_data` to have mean 0 and standard deviation 1.0")
  }

  # return stuff
  Return = list( "X_gtp"=X_gtp, "X_itp"=X_itp, "covariate_names"=covariate_names )
  return( Return )
}

