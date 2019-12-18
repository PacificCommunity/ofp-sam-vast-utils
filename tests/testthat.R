library(testthat)
library(ndd.vast.utils)

# set path to where the data is
	testdata.dir = system.file("testdata", package="ndd.vast.utils")

test_check("ndd.vast.utils")
