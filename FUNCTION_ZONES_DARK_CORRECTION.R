################################################################################################
################################################################################################
# Function for define zones for new dark correction in OMZs and STGs (+ Baltic and Black seas) #
################################################################################################
###########################  Raphaëlle Sauzède  - August 2020 ##################################
################################################################################################
################################################################################################

library(ncdf4)
library(OceanView)

# path to WOA 2018data
# Data downloaded from https://www.nodc.noaa.gov/cgi-bin/OC5/woa18/woa18oxnu.pl the 14/08/2020 (1° grid)
data_woa_doxy <- "/home/sauzede/Documents/R/DATA/WOA/DOXY/woa18_all_o00_01.nc"

# path to NASA entire composite Chl climatological data
# Data downloaded from https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A20021852020182.L3m_CU_CHL_chlor_a_4km.nc (1/24° grid) 
# the 17/08/2020
data_nasa_chl <- "/home/sauzede/Documents/R/DATA/OC/NASA/CHL/COMPOSITE_CLIMATO/A20021852020182.L3m_CU_CHL_chlor_a_4km.nc"

#Open the DOXY NetCDF file and get lon/lat/depth
#o_an : Objectively analyzed mean fields for mole_concentration_of_dissolved_molecular_oxygen_in_sea_water at standard depth levels
nc_doxy <- nc_open(data_woa_doxy)
doxy_data <- ncvar_get(nc_doxy, "o_an")
lon_doxy_data <- ncvar_get(nc_doxy, "lon")
lat_doxy_data <- ncvar_get(nc_doxy, "lat")
depth_data <- ncvar_get(nc_doxy, "depth")
# Choose the depth 400 (for OMZ)
doxy_data <- doxy_data[,,which(depth_data==400)]
nc_close(nc_doxy)

#Open the Chl NetCDF file and get lon/lat/depth
nc_chl <- nc_open(data_nasa_chl)
chl_data <- ncvar_get(nc_chl, "chlor_a")
lon_chl_data <- ncvar_get(nc_chl, "lon")
lat_chl_data <- ncvar_get(nc_chl, "lat")
chl_data <- chl_data[order(lon_chl_data), order(lat_chl_data)]
lon_chl_data <- lon_chl_data[order(lon_chl_data)]
lat_chl_data <- lat_chl_data[order(lat_chl_data)]
nc_close(nc_chl)

# Get a coarser version of Chl sat in order to have more homogeneous oligotrophic zones (in order to smooth because of some pixels at high resolution)
# Use the function remap to get chl in a 0.5° grid
lon_chl_data2 <- seq(-179.5,179.5,0.5)
lat_chl_data2 <- seq(-89.5,89.5,0.5)
# Remap the Chl on the new grid
chl_data_coarser <- remap(chl_data, x = lon_chl_data, y = lat_chl_data, xto= lon_chl_data2, yto = lat_chl_data2)$var

#######################################################################################################
############################    MAP CODE  #############################################################
## Map doxy data from WOA at 400 m  depth with the 70 umol kg-1 delimitation
# library(fields)
# library(maps)
# png("/home/sauzede/Documents/R/PLOT/SOCA_BBP/SOCA_2020/SOCA_CHL_RENOSH/WOA_DOXY_400.png", res=300, width=400*7, height= 400*4)
# image.plot(x = lon_doxy_data, y = lat_doxy_data, doxy_data, xlab = "Longitude", ylab = "Latitude")
# contour(x = lon_doxy_data, y = lat_doxy_data, doxy_data, levels = c(70), add=T)
# map(add=T, col="grey10", fill = T, lty=0)
# graphics.off()

# # Map chl data from NASA climatology (entire composite) with the 0.1 mg m-3 delimitation
# library(fields)
# library(maps)
# # with high resolution data
# # png("/home/sauzede/Documents/R/PLOT/SOCA_BBP/SOCA_2020/SOCA_CHL_RENOSH/NASA_CHL_COARSER_0.08.png", res=300, width=400*7, height= 400*4)
# # image.plot(x = lon_chl_data, y = lat_chl_data, log10(chl_data), xlab = "Longitude", ylab = "Latitude", axis.args = list(at=log10(c(0.01,0.1,1,10,50)), labels = c(0.01,0.1,1,10,50)))
# # contour(x = lon_chl_data, y = lat_chl_data, chl_data, levels = c(0.08), add=T)
# # with low resolution data
# png("/home/sauzede/Documents/R/PLOT/SOCA_BBP/SOCA_2020/SOCA_CHL_RENOSH/NASA_CHL_COARSER_0.08.png", res=300, width=400*7, height= 400*4)
# image.plot(x = lon_chl_data2, y = lat_chl_data2, log10(chl_data_coarser), xlab = "Longitude", ylab = "Latitude", axis.args = list(at=log10(c(0.01,0.1,1,10,50)), labels = c(0.01,0.1,1,10,50)))
# contour(x = lon_chl_data2, y = lat_chl_data2, chl_data_coarser, levels = c(0.08), add=T)
# 
# # IME zone near Hawai
# rect(xleft = -175, xright = -150, ybottom = 18, ytop = 28)
# # IME zone near Tahiti and Marquesas Islands
# rect(xleft = -155, xright = -130, ybottom = -25, ytop = -10)
# # IME zone near Samoa
# rect(xleft = -175, xright = -170, ybottom = -15, ytop = -12)
# # IME zone near Wallis Futuna, Vanuatu, Nouvelle calédonie etc.
# rect(xleft = 125, xright = 180, ybottom = 5, ytop = 12)
# 
# map(add=T, col="grey10", fill = T, lty=0)
# graphics.off()

##########################################################################################################
##########################################################################################################

# The threshold for doxy is set initially to 75 µmol kg-1
# The threshold for Chl is set initially to 0.08 mg m-3
Zone <- function(lat, lon, threshold_doxy = 75, threshold_chl = 0.08){
  
  zone <- "OTHER"

  # Get the closest pixel for doxy_data
  iii_lon_doxy <- which.min(abs(lon_doxy_data-lon))
  iii_lat_doxy <- which.min(abs(lat_doxy_data-lat))
  # Compute the mean with closest pixels because of NaN values
  doxy <- mean(doxy_data[seq(max(1,iii_lon_doxy-2), min(iii_lon_doxy+2,360)),seq(max(1,iii_lat_doxy-2), min(iii_lat_doxy+2,180))], na.rm=T)
  if(!is.na(doxy) & doxy<=threshold_doxy){zone <- "OMZ"}
  
  # Get the closest pixel for chl_data because chl data is almost non empty 
  iii_lon_chl <- which.min(abs(lon_chl_data2-lon))
  iii_lat_chl <- which.min(abs(lat_chl_data2-lat))
  chlor_a <- chl_data_coarser[iii_lon_chl,iii_lat_chl]
  if(!is.na(chlor_a) & chlor_a<=threshold_chl){zone <- "SUBTROPICAL_GYRE"}
  
  # When there are island mass effects (i.e. the higher chl values closer to islands) chl is < threshold but we want to correct these data!
  # So we add manually some areas to remove
  # IME zone near Hawai
  if(lon > -175 & lon < -150 & lat > 18 & lat < 28){zone <- "SUBTROPICAL_GYRE"}
  # IME zone near Tahiti and Marquesas Islands
  if(lon > -155 & lon < -130 & lat > -25 &  lat < -10){zone <- "SUBTROPICAL_GYRE"}
  # IME zone near Samoa
  if(lon > -175  & lon < -170 & lat > -15 & lat < -12){zone <- "SUBTROPICAL_GYRE"}
  # IME zone near Wallis Futuna, Vanuatu, Nouvelle calédonie etc.
  if(lon > 125 & lon < 180 & lat > 5 & lat < 12){zone <- "SUBTROPICAL_GYRE"}
  
  # Baltic sea to remove
  lon_min_BaltS <- 9
  lon_max_BaltS <- 23
  lat_min_BaltS <- 52
  lat_max_BaltS <- 65
  if(lon>lon_min_BaltS & lon<lon_max_BaltS & lat>lat_min_BaltS & lat<lat_max_BaltS){zone <- "BALTIC_SEA"}
  
  # Black sea to remove
  lon_min_BlackS <- 25
  lon_max_BlackS <- 41
  lat_min_BlackS <- 40
  lat_max_BlackS <- 47
  if(lon>lon_min_BlackS & lon<lon_max_BlackS & lat>lat_min_BlackS & lat<lat_max_BlackS){zone <- "BLACK_SEA"}
  
  return(zone)
}
