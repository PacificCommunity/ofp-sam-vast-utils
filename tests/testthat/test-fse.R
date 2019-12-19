test_that("fine scale, strata, and enviro is working",
	{
		skip_on_travis()
  		# Bring in test data
  			data(Data_Geostat)
  			data(enviro)
  			data(skj.alt2019.shp)

  		# set path to where the data is
			testdata.dir = system.file("testdata", package="ndd.vast.utils")

  		# Bring in solution
  			load( file.path(testdata.dir,"/fse.mgc.RData") )
  			load( file.path(testdata.dir,"/fse.par.RData") )
  			load( file.path(testdata.dir,"/fse.index_cyl.RData") )

  		# fit model
  			vast_output = fit.vast(Data_Geostat,RunDir=getwd(),Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=TRUE,input.grid.res=1,knot_method = "grid",n_x=100,Version="VAST_v8_3_0",Method="Mesh",strata.sp=skj.alt2019.shp,enviro=enviro)
  			Par1 = fse.par
  			Par2 = vast_output$Opt$par

  		expect_equal( as.vector(Par1), as.vector(Par2), tolerance=1e-3 )
		expect_equal( as.vector(fse.index_cyl), as.vector(vast_output$Report$Index_cyl), tolerance=1e-2 )
		expect_lt(fse.mgc,expected=1e-4)
	})
