
#' Function to convert points stored as eastings and northings to lat-lon
#' 
#' @param E Vector of eastings
#' @param N Vector of northings
#' @param crs.en Character string of the crs for the E-N projection
#' @param crs.ll Character string of the crs for the current lat-lon projections
#' @return Matrix of two columns: Lon and Lat
#' @importFrom sp coordinates
#' @importFrom sp proj4string
#' @importFrom sp CRS
#' @importFrom sp spTransform
#' @keywords internal

Convert_EN_to_LL_Fn.ndd = function(E, N, crs.en = "+proj=tpeqd +lat_1=0 +lon_1=155 +lat_2=0 +lon_2=209 +datum=WGS84 +ellps=WGS84 +units=km +no_defs",crs.ll = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# modified from Jim Thorson's FishStatUtils::Convert_LL_to_EastNorth_Fn.R
# pass x and y locations as Lon and Lat
# crs default is two point equidistant # OPTIONAL: +a=6005 appears to be the distance in km between the two points
{
  # Transform
      dstart=data.frame(E=E, N=N) # that's the object
      sp::coordinates(dstart) = c("E", "N")
      sp::proj4string(dstart) = sp::CRS(crs.en) # that's the lat long projection
      CRS.new = sp::CRS(crs.ll) # that's the eastings and northings projection
      dstart.t = sp::spTransform(dstart, CRS.new) # here's where you transform

  # Clean up
      dstart.t = cbind( "Lon"=dstart.t@coords[,"E"], "Lat"=dstart.t@coords[,'N'])
      dstart.t[,1] = ifelse(dstart.t[,1]<0,dstart.t[,1]+360,dstart.t)
      attr(dstart.t,"zone") = crs.ll

  # Return results
    return( dstart.t )
}

#' Function to convert points stored as lat-lon to eastings and northings
#' 
#' @param Lon Vector of longitudes
#' @param Lat Vector of latitudes
#' @param crs.en Character string of the crs for the E-N projection
#' @param crs.ll Character string of the crs for the current lat-lon projections
#' @return Matrix of two columns: E_km and N_km
#' @importFrom sp coordinates
#' @importFrom sp proj4string
#' @importFrom sp CRS
#' @importFrom sp spTransform
#' @keywords internal

Convert_LL_to_EastNorth_Fn.ndd = function(Lon, Lat, crs.en = "+proj=tpeqd +lat_1=0 +lon_1=155 +lat_2=0 +lon_2=209 +datum=WGS84 +ellps=WGS84 +units=km +no_defs",crs.ll = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# modified from Jim Thorson's FishStatUtils::Convert_LL_to_EastNorth_Fn.R
# pass x and y locations as Lon and Lat
# crs default is two point equidistant # OPTIONAL: +a=6005 appears to be the distance in km between the two points
{

  # Transform
      dstart=data.frame(lon=Lon, lat=Lat) # that's the object
      sp::coordinates(dstart) = c("lon", "lat")
      sp::proj4string(dstart) = sp::CRS(crs.ll) # that's the lat long projection
      CRS.new = sp::CRS(crs.en) # that's the eastings and northings projection
      dstart.t = sp::spTransform(dstart, CRS.new) # here's where you transform

  # Clean up
      dstart.t = cbind( "E_km"=dstart.t@coords[,"lon"], "N_km"=dstart.t@coords[,'lat'])
      attr(dstart.t,"zone") = crs.en

  # Return results
    return( dstart.t )
}

#' # code to create "turbo" colormap
#' All credit to those listed below; this is just a reposting of their work to utilize the colorscale within the package
#' Better, shorter, more idiomatic, vectorized versions of the interpolate and turbo functions 
#' thanks to user 'onesandzeroes'
#' All credit to those listed below; this is just a reposting of their work to utilize the colorscale within the package
#' https://ai.googleblog.com/2019/08/turbo-improved-rainbow-colormap-for.html
#' https://gist.github.com/mikhailov-work/ee72ba4191942acecc03fe6da94fc73f
#' https://gist.github.com/jlmelville/be981e2f36485d8ef9616aef60fd52ab
#' 
#' @param colormap ...
#' @param x ...
#' @keywords internal
#' @noRd
 
interpolate_vec <- function(colormap, x) {
  x <- pmax(0.0, pmin(1.0, x))
  a <- floor(x * 255.0)
  b <- pmin(255, a + 1)
  f <- x * 255.0 - a
  a <- a + 1
  b <- b + 1
  
  colormap[a, ] + (colormap[b, ] - colormap[a, ]) * f
}

#' # code to create "turbo" colormap
#' All credit to those listed below; this is just a reposting of their work to utilize the colorscale within the package
#' https://ai.googleblog.com/2019/08/turbo-improved-rainbow-colormap-for.html
#' https://gist.github.com/mikhailov-work/ee72ba4191942acecc03fe6da94fc73f
#' https://gist.github.com/jlmelville/be981e2f36485d8ef9616aef60fd52ab
#' 
#' @param n Number of colors from 'turbo' palette to generate
#' @param start Where on scale to start [0,1]
#' @param end Where on scale to end [0,1]
#' @keywords internal
#' @noRd


turbo_vec <- function(n, start = 0, end = 1) {
      turbo_colormap_data <- matrix(
      c(
        c(0.18995, 0.07176, 0.23217),
        c(0.19483, 0.08339, 0.26149),
        c(0.19956, 0.09498, 0.29024),
        c(0.20415, 0.10652, 0.31844),
        c(0.20860, 0.11802, 0.34607),
        c(0.21291, 0.12947, 0.37314),
        c(0.21708, 0.14087, 0.39964),
        c(0.22111, 0.15223, 0.42558),
        c(0.22500, 0.16354, 0.45096),
        c(0.22875, 0.17481, 0.47578),
        c(0.23236, 0.18603, 0.50004),
        c(0.23582, 0.19720, 0.52373),
        c(0.23915, 0.20833, 0.54686),
        c(0.24234, 0.21941, 0.56942),
        c(0.24539, 0.23044, 0.59142),
        c(0.24830, 0.24143, 0.61286),
        c(0.25107, 0.25237, 0.63374),
        c(0.25369, 0.26327, 0.65406),
        c(0.25618, 0.27412, 0.67381),
        c(0.25853, 0.28492, 0.69300),
        c(0.26074, 0.29568, 0.71162),
        c(0.26280, 0.30639, 0.72968),
        c(0.26473, 0.31706, 0.74718),
        c(0.26652, 0.32768, 0.76412),
        c(0.26816, 0.33825, 0.78050),
        c(0.26967, 0.34878, 0.79631),
        c(0.27103, 0.35926, 0.81156),
        c(0.27226, 0.36970, 0.82624),
        c(0.27334, 0.38008, 0.84037),
        c(0.27429, 0.39043, 0.85393),
        c(0.27509, 0.40072, 0.86692),
        c(0.27576, 0.41097, 0.87936),
        c(0.27628, 0.42118, 0.89123),
        c(0.27667, 0.43134, 0.90254),
        c(0.27691, 0.44145, 0.91328),
        c(0.27701, 0.45152, 0.92347),
        c(0.27698, 0.46153, 0.93309),
        c(0.27680, 0.47151, 0.94214),
        c(0.27648, 0.48144, 0.95064),
        c(0.27603, 0.49132, 0.95857),
        c(0.27543, 0.50115, 0.96594),
        c(0.27469, 0.51094, 0.97275),
        c(0.27381, 0.52069, 0.97899),
        c(0.27273, 0.53040, 0.98461),
        c(0.27106, 0.54015, 0.98930),
        c(0.26878, 0.54995, 0.99303),
        c(0.26592, 0.55979, 0.99583),
        c(0.26252, 0.56967, 0.99773),
        c(0.25862, 0.57958, 0.99876),
        c(0.25425, 0.58950, 0.99896),
        c(0.24946, 0.59943, 0.99835),
        c(0.24427, 0.60937, 0.99697),
        c(0.23874, 0.61931, 0.99485),
        c(0.23288, 0.62923, 0.99202),
        c(0.22676, 0.63913, 0.98851),
        c(0.22039, 0.64901, 0.98436),
        c(0.21382, 0.65886, 0.97959),
        c(0.20708, 0.66866, 0.97423),
        c(0.20021, 0.67842, 0.96833),
        c(0.19326, 0.68812, 0.96190),
        c(0.18625, 0.69775, 0.95498),
        c(0.17923, 0.70732, 0.94761),
        c(0.17223, 0.71680, 0.93981),
        c(0.16529, 0.72620, 0.93161),
        c(0.15844, 0.73551, 0.92305),
        c(0.15173, 0.74472, 0.91416),
        c(0.14519, 0.75381, 0.90496),
        c(0.13886, 0.76279, 0.89550),
        c(0.13278, 0.77165, 0.88580),
        c(0.12698, 0.78037, 0.87590),
        c(0.12151, 0.78896, 0.86581),
        c(0.11639, 0.79740, 0.85559),
        c(0.11167, 0.80569, 0.84525),
        c(0.10738, 0.81381, 0.83484),
        c(0.10357, 0.82177, 0.82437),
        c(0.10026, 0.82955, 0.81389),
        c(0.09750, 0.83714, 0.80342),
        c(0.09532, 0.84455, 0.79299),
        c(0.09377, 0.85175, 0.78264),
        c(0.09287, 0.85875, 0.77240),
        c(0.09267, 0.86554, 0.76230),
        c(0.09320, 0.87211, 0.75237),
        c(0.09451, 0.87844, 0.74265),
        c(0.09662, 0.88454, 0.73316),
        c(0.09958, 0.89040, 0.72393),
        c(0.10342, 0.89600, 0.71500),
        c(0.10815, 0.90142, 0.70599),
        c(0.11374, 0.90673, 0.69651),
        c(0.12014, 0.91193, 0.68660),
        c(0.12733, 0.91701, 0.67627),
        c(0.13526, 0.92197, 0.66556),
        c(0.14391, 0.92680, 0.65448),
        c(0.15323, 0.93151, 0.64308),
        c(0.16319, 0.93609, 0.63137),
        c(0.17377, 0.94053, 0.61938),
        c(0.18491, 0.94484, 0.60713),
        c(0.19659, 0.94901, 0.59466),
        c(0.20877, 0.95304, 0.58199),
        c(0.22142, 0.95692, 0.56914),
        c(0.23449, 0.96065, 0.55614),
        c(0.24797, 0.96423, 0.54303),
        c(0.26180, 0.96765, 0.52981),
        c(0.27597, 0.97092, 0.51653),
        c(0.29042, 0.97403, 0.50321),
        c(0.30513, 0.97697, 0.48987),
        c(0.32006, 0.97974, 0.47654),
        c(0.33517, 0.98234, 0.46325),
        c(0.35043, 0.98477, 0.45002),
        c(0.36581, 0.98702, 0.43688),
        c(0.38127, 0.98909, 0.42386),
        c(0.39678, 0.99098, 0.41098),
        c(0.41229, 0.99268, 0.39826),
        c(0.42778, 0.99419, 0.38575),
        c(0.44321, 0.99551, 0.37345),
        c(0.45854, 0.99663, 0.36140),
        c(0.47375, 0.99755, 0.34963),
        c(0.48879, 0.99828, 0.33816),
        c(0.50362, 0.99879, 0.32701),
        c(0.51822, 0.99910, 0.31622),
        c(0.53255, 0.99919, 0.30581),
        c(0.54658, 0.99907, 0.29581),
        c(0.56026, 0.99873, 0.28623),
        c(0.57357, 0.99817, 0.27712),
        c(0.58646, 0.99739, 0.26849),
        c(0.59891, 0.99638, 0.26038),
        c(0.61088, 0.99514, 0.25280),
        c(0.62233, 0.99366, 0.24579),
        c(0.63323, 0.99195, 0.23937),
        c(0.64362, 0.98999, 0.23356),
        c(0.65394, 0.98775, 0.22835),
        c(0.66428, 0.98524, 0.22370),
        c(0.67462, 0.98246, 0.21960),
        c(0.68494, 0.97941, 0.21602),
        c(0.69525, 0.97610, 0.21294),
        c(0.70553, 0.97255, 0.21032),
        c(0.71577, 0.96875, 0.20815),
        c(0.72596, 0.96470, 0.20640),
        c(0.73610, 0.96043, 0.20504),
        c(0.74617, 0.95593, 0.20406),
        c(0.75617, 0.95121, 0.20343),
        c(0.76608, 0.94627, 0.20311),
        c(0.77591, 0.94113, 0.20310),
        c(0.78563, 0.93579, 0.20336),
        c(0.79524, 0.93025, 0.20386),
        c(0.80473, 0.92452, 0.20459),
        c(0.81410, 0.91861, 0.20552),
        c(0.82333, 0.91253, 0.20663),
        c(0.83241, 0.90627, 0.20788),
        c(0.84133, 0.89986, 0.20926),
        c(0.85010, 0.89328, 0.21074),
        c(0.85868, 0.88655, 0.21230),
        c(0.86709, 0.87968, 0.21391),
        c(0.87530, 0.87267, 0.21555),
        c(0.88331, 0.86553, 0.21719),
        c(0.89112, 0.85826, 0.21880),
        c(0.89870, 0.85087, 0.22038),
        c(0.90605, 0.84337, 0.22188),
        c(0.91317, 0.83576, 0.22328),
        c(0.92004, 0.82806, 0.22456),
        c(0.92666, 0.82025, 0.22570),
        c(0.93301, 0.81236, 0.22667),
        c(0.93909, 0.80439, 0.22744),
        c(0.94489, 0.79634, 0.22800),
        c(0.95039, 0.78823, 0.22831),
        c(0.95560, 0.78005, 0.22836),
        c(0.96049, 0.77181, 0.22811),
        c(0.96507, 0.76352, 0.22754),
        c(0.96931, 0.75519, 0.22663),
        c(0.97323, 0.74682, 0.22536),
        c(0.97679, 0.73842, 0.22369),
        c(0.98000, 0.73000, 0.22161),
        c(0.98289, 0.72140, 0.21918),
        c(0.98549, 0.71250, 0.21650),
        c(0.98781, 0.70330, 0.21358),
        c(0.98986, 0.69382, 0.21043),
        c(0.99163, 0.68408, 0.20706),
        c(0.99314, 0.67408, 0.20348),
        c(0.99438, 0.66386, 0.19971),
        c(0.99535, 0.65341, 0.19577),
        c(0.99607, 0.64277, 0.19165),
        c(0.99654, 0.63193, 0.18738),
        c(0.99675, 0.62093, 0.18297),
        c(0.99672, 0.60977, 0.17842),
        c(0.99644, 0.59846, 0.17376),
        c(0.99593, 0.58703, 0.16899),
        c(0.99517, 0.57549, 0.16412),
        c(0.99419, 0.56386, 0.15918),
        c(0.99297, 0.55214, 0.15417),
        c(0.99153, 0.54036, 0.14910),
        c(0.98987, 0.52854, 0.14398),
        c(0.98799, 0.51667, 0.13883),
        c(0.98590, 0.50479, 0.13367),
        c(0.98360, 0.49291, 0.12849),
        c(0.98108, 0.48104, 0.12332),
        c(0.97837, 0.46920, 0.11817),
        c(0.97545, 0.45740, 0.11305),
        c(0.97234, 0.44565, 0.10797),
        c(0.96904, 0.43399, 0.10294),
        c(0.96555, 0.42241, 0.09798),
        c(0.96187, 0.41093, 0.09310),
        c(0.95801, 0.39958, 0.08831),
        c(0.95398, 0.38836, 0.08362),
        c(0.94977, 0.37729, 0.07905),
        c(0.94538, 0.36638, 0.07461),
        c(0.94084, 0.35566, 0.07031),
        c(0.93612, 0.34513, 0.06616),
        c(0.93125, 0.33482, 0.06218),
        c(0.92623, 0.32473, 0.05837),
        c(0.92105, 0.31489, 0.05475),
        c(0.91572, 0.30530, 0.05134),
        c(0.91024, 0.29599, 0.04814),
        c(0.90463, 0.28696, 0.04516),
        c(0.89888, 0.27824, 0.04243),
        c(0.89298, 0.26981, 0.03993),
        c(0.88691, 0.26152, 0.03753),
        c(0.88066, 0.25334, 0.03521),
        c(0.87422, 0.24526, 0.03297),
        c(0.86760, 0.23730, 0.03082),
        c(0.86079, 0.22945, 0.02875),
        c(0.85380, 0.22170, 0.02677),
        c(0.84662, 0.21407, 0.02487),
        c(0.83926, 0.20654, 0.02305),
        c(0.83172, 0.19912, 0.02131),
        c(0.82399, 0.19182, 0.01966),
        c(0.81608, 0.18462, 0.01809),
        c(0.80799, 0.17753, 0.01660),
        c(0.79971, 0.17055, 0.01520),
        c(0.79125, 0.16368, 0.01387),
        c(0.78260, 0.15693, 0.01264),
        c(0.77377, 0.15028, 0.01148),
        c(0.76476, 0.14374, 0.01041),
        c(0.75556, 0.13731, 0.00942),
        c(0.74617, 0.13098, 0.00851),
        c(0.73661, 0.12477, 0.00769),
        c(0.72686, 0.11867, 0.00695),
        c(0.71692, 0.11268, 0.00629),
        c(0.70680, 0.10680, 0.00571),
        c(0.69650, 0.10102, 0.00522),
        c(0.68602, 0.09536, 0.00481),
        c(0.67535, 0.08980, 0.00449),
        c(0.66449, 0.08436, 0.00424),
        c(0.65345, 0.07902, 0.00408),
        c(0.64223, 0.07380, 0.00401),
        c(0.63082, 0.06868, 0.00401),
        c(0.61923, 0.06367, 0.00410),
        c(0.60746, 0.05878, 0.00427),
        c(0.59550, 0.05399, 0.00453),
        c(0.58336, 0.04931, 0.00486),
        c(0.57103, 0.04474, 0.00529),
        c(0.55852, 0.04028, 0.00579),
        c(0.54583, 0.03593, 0.00638),
        c(0.53295, 0.03169, 0.00705),
        c(0.51989, 0.02756, 0.00780),
        c(0.50664, 0.02354, 0.00863),
        c(0.49321, 0.01963, 0.00955),
        c(0.47960, 0.01583, 0.01055)
      ),
      ncol = 3,
      byrow = TRUE
    )
  xs <- seq.int(from = start, to = end, length.out = n)
  rgb(interpolate_vec(turbo_colormap_data, xs))
}
