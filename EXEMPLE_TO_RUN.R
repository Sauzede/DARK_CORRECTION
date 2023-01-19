###############################################################################################
###############################################################################################
#################  Code to run as example the Dark correction function ########################
######  And the function to find zones where this correction has to be applied ################
###############################################################################################
#########################  Raphaëlle Sauzède - February 2021 ##################################
###############################################################################################
###############################################################################################
###############################################################################################

# Library to load
library(ncdf4)
library(maps)
library(oce)
library(fields)
library(viridis)
library(caTools)
library(stringr)


# Define working path data
path_data_home <- "/home/sauzede/Documents/R/DATA/"
path_data_argo <- "/home/sauzede/ARGO_DATA/"

setwd(dir = "/home/sauzede/Documents/R/CODE/FONCTIONS/CODE_DARK_CORRECTION_CLEAN/")
path_function <- "FUNCTIONS/"
path_plot <- "PLOT/"

# Source the functions
source(paste(path_function,"FUNCTION_ZONES_DARK_CORRECTION.R",sep=""))
source(paste(path_function,"FUNCTION_DARK_CORRECTION.R",sep=""))
source(paste(path_function,"FUNCTION_INTERP_NAN.R",sep=""))
source(paste(path_function,"FUNCTION_RUNNING_FILTER.R",sep=""))

date_today = "20210218"

#################################################################################################################
################################     SELECT DATA FILES WITH CHL   ###############################################
#################################################################################################################

FILE_INDEX <- paste0(path_data_home,"/index_bio_",date_today,".txt")

if(!file.exists(FILE_INDEX)){
  # download file with information on Real Time/ Adjusted and Delayed mode parameters
  download.file("ftp://ftp.ifremer.fr/ifremer/argo/argo_bio-profile_index.txt",
                FILE_INDEX, quiet = FALSE, mode = "w",
                cacheOK = TRUE)}


#read the Index file
INDEX_BIO <- read.table(FILE_INDEX,skip=8,h=T,sep=",")

#Total number of profiles
length(INDEX_BIO[,1])


files <- as.character(INDEX_BIO$file) #pathways of profiles
ident <- strsplit(files,"/") #split DAC/WMO and Netcdf files
ident <- matrix(unlist(ident), ncol=4, byrow=TRUE) # Put it in a matrix
dac <- ident[,1] # Derive dac
wod <- ident[,2] # Derive wmo
ncdf <- ident[,4] # Derive ncdf file
variables<-as.character(INDEX_BIO$parameters) #contains parameters measured by the float
mode <- as.character(INDEX_BIO$parameter_data_mode) #contains the mode for each parameter (R for real time, A for adjusted and D for Delayed mode)
# Get lon/lat/date
lat <- INDEX_BIO$latitude
lon <- INDEX_BIO$longitude
#ref_date <- "19500101000000"
INDEX_BIO$date <- as.character(as.numeric(INDEX_BIO$date)) 
date <- as.POSIXct(INDEX_BIO$date, format="%Y%m%d%H%M%S")


variables <- strsplit(variables," ") #separate the different variables
## selector for the selection of profiles with BGC PARAM wanted measured (here CHLA)
selector <- NA
for (i in 1:length(files)) {
  print(paste(i,"/",length(files)))
  # verify that CHLA is measured
  if ("CHLA" %in% variables[[i]]){
    selector[i]<-1
  }
  else{
    selector[i] <- 0
  }
  # verify that it is not a descendant profile
  n_cycle <- unlist(strsplit(unlist(strsplit(unlist(strsplit(files[[i]],"/"))[4],".nc")),"_"))[2]
  if(length(unlist(strsplit(n_cycle,"")))>3){
    selector[i] <- 0
  }
}

print(paste(length(selector[selector==1]), "files with CHLA measured on a total of", length(selector), "files", sep = " "))

# Get the index and metadata only for the files with CHLA measured
INDEX_BIO_SELECT <- INDEX_BIO[selector==1,]

# Get the paths to CHLA files
chemin_files <- as.character(INDEX_BIO_SELECT$file)
chemin_files <- chemin_files[!is.na(chemin_files)]

# Get zones (application of dark correction or not) from lon/lat of profiles
# Get lon/lat
lat_select <- INDEX_BIO_SELECT$latitude
lon_select <- INDEX_BIO_SELECT$longitude
# remove NA position
lat_select <- lat_select[!is.na(lon_select)]
chemin_files <- chemin_files[!is.na(lon_select)]
lon_select <- lon_select[!is.na(lon_select)]
  
# Apply function to determine the zones with thresholds
# Threshold for DOXY is initially set to 75 µmol kg-1
# Threshold for Chl to define subtropical gyres is initially set to 0.08 mg m-3
zones_select <- NULL
for(i in 1:length(lon_select)){
  zones_select <- c(zones_select,Zone(lat = lat_select[i], lon = lon_select[i]))
  print(i)
}

unique(zones_select)

####### Geolocation map of CHLA profiles with associated zones
# Load chl and doxy data to verify the identification of the different zones (with contour)
data_woa_doxy <- paste0(path_data_home, "WOA/DOXY/woa18_all_o00_01.nc")
data_nasa_chl <- paste0(path_data_home,"OC/NASA/CHL/COMPOSITE_CLIMATO/A20021852020182.L3m_CU_CHL_chlor_a_4km.nc")

#Open the DOXY NetCDF file and get lon/lat/depth and DOXY at 400 m depth
#o_an : Objectively analyzed mean fields for mole_concentration_of_dissolved_molecular_oxygen_in_sea_water at standard depth levels
nc_doxy <- nc_open(data_woa_doxy) # Open the NetCDF file
doxy_data <- ncvar_get(nc_doxy, "o_an") # Get the DOXY data
lon_doxy_data <- ncvar_get(nc_doxy, "lon") # Get longitude
lat_doxy_data <- ncvar_get(nc_doxy, "lat") # Get latitude
depth_data <- ncvar_get(nc_doxy, "depth") # Get depth
doxy_data <- doxy_data[,,which(depth_data==400)] # Choose the depth 400 (for OMZ)
nc_close(nc_doxy) # Close the NetCDF file 

#Open the Chl NetCDF file and get lon/lat/depth and Chl
nc_chl <- nc_open(data_nasa_chl) # Open the NetCDF file
chl_data <- ncvar_get(nc_chl, "chlor_a") # Get chlorophyll data
lon_chl_data <- ncvar_get(nc_chl, "lon") # Get longitude
lat_chl_data <- ncvar_get(nc_chl, "lat") # Get latitude
chl_data <- chl_data[order(lon_chl_data), order(lat_chl_data)] # Order  matrixdata with lon/lat
lon_chl_data <- lon_chl_data[order(lon_chl_data)] # Order vector data with lon
lat_chl_data <- lat_chl_data[order(lat_chl_data)] # Order vector data with lat
nc_close(nc_chl) # Close the NetCDF file


# We use a coarser resolution for Chla to define more homogeneous zones
# Use the function remap to get chl in a 0.5° grid
lon_chl_data2 <- seq(-179.5,179.5,0.5)
lat_chl_data2 <- seq(-89.5,89.5,0.5)
# Remap the data using the new grid
chl_data_coarser <- remap(chl_data, x = lon_chl_data, y = lat_chl_data, xto= lon_chl_data2, yto = lat_chl_data2)$var

png(paste0(path_plot,"profiles_with_zones.png"),res=300, width=400*7, height= 400*4)
# Plot the chl data and the contour of doxy (dashed line) and chl threshold
image.plot(x = lon_chl_data, y = lat_chl_data, log10(chl_data), xlab = "Longitude", ylab = "Latitude", axis.args = list(at=log10(c(0.01,0.1,1,10,50)), labels = c(0.01,0.1,1,10,50)))
contour(x = lon_doxy_data, y = lat_doxy_data, doxy_data, levels = c(75), add=T, lty =2)
contour(x = lon_chl_data2, y = lat_chl_data2, chl_data_coarser, levels = c(0.08), add=T)
contour(x = lon_chl_data2, y = lat_chl_data2, chl_data_coarser, levels = c(0.08), add=T)

# Add profiles for different zones
points(lon_select[zones_select == "OMZ"], lat_select[zones_select == "OMZ"], 
       pch=19, col="orange",
       ylim=c(-80,80), xlim=c(-170,170),
       xlab = "Longitude", ylab= "Latitude", cex=.5)
points(lon_select[zones_select =="SUBTROPICAL_GYRE"], lat_select[zones_select == "SUBTROPICAL_GYRE"], 
       pch=19, col="turquoise4")
points(lon_select[zones_select =="BALTIC_SEA"], lat_select[zones_select == "BALTIC_SEA"], 
       pch=19, col="chartreuse3")
points(lon_select[zones_select =="BLACK_SEA"], lat_select[zones_select == "BLACK_SEA"], 
       pch=19, col="violetred4")
points(lon_select[zones_select =="OTHER"], lat_select[zones_select == "OTHER"], 
       pch=19, col="grey", cex=.3)
map(add=T, col="grey10", fill = T, lty=0)
legend("bottomright", legend = c("OMZ", "Oligotrophic subtropical gyres", "Black Sea", "Baltic Sea"), 
       col = c("orange", "turquoise4", "violetred4", "chartreuse3"), box.lty = 0, pch = 19, ncol = 2, bg = "white", cex=.7)
graphics.off()

# Plot the BGC-Argo profiles that are in the zones to correct from dark specifically
png(paste0(path_plot,"profiles_zones_dark_correction.png"), res=300, width=400*7, height= 400*4)
plot(lon_select[zones_select == "OMZ"], lat_select[zones_select == "OMZ"], 
     pch=19, col="orange",
     ylim=c(-80,80), xlim=c(-170,170),
     xlab = "Longitude", ylab= "Latitude", cex=.5)
points(lon_select[zones_select =="SUBTROPICAL_GYRE"], lat_select[zones_select == "SUBTROPICAL_GYRE"], 
       pch=19, col="turquoise4")
points(lon_select[zones_select =="BALTIC_SEA"], lat_select[zones_select == "BALTIC_SEA"], 
       pch=19, col="chartreuse3")
points(lon_select[zones_select =="BLACK_SEA"], lat_select[zones_select == "BLACK_SEA"], 
       pch=19, col="violetred4")
points(lon_select[zones_select =="OTHER"], lat_select[zones_select == "OTHER"], 
       pch=19, col="grey", cex=.3)
map(add=T, col="grey10", fill = T, lty=0)
legend("bottomright", legend = c("OMZ", "Oligotrophic subtropical gyres", "Black Sea", "Baltic Sea"), 
       col = c("orange", "turquoise4", "violetred4", "chartreuse3"), box.lty = 0, pch = 19, ncol = 2, bg = "white", cex=.7)
graphics.off()

# Get path to files for dark correction --> OMZ zones + STG + Baltic Sea + Black Sea
chemin_files_to_correct <- chemin_files[zones_select %in% c("OMZ", "SUBTROPICAL_GYRE", "BALTIC_SEA", "BLACK_SEA")]

print(paste(length(chemin_files_to_correct), "CHLA profiles to apply the specific dark correction on a total of", length(selector), "files", sep = " "))

# Apply the new dark correction for these and plot it in a pdf file
pdf(paste(path_plot,"profiles_CHL_correction_Wojtasiewicz.pdf",sep=""))
par(mfrow=c(2,2))
par(las=1)
# Loop on all profiles to correct
for(c in 1:length(chemin_files_to_correct)){
  
  # Get dac and wmo from the path
  ident <- strsplit(chemin_files_to_correct[c],"/") #split DAC/WMO and Netcdf files
  ident <- matrix(unlist(ident), ncol=4, byrow=TRUE) # Get the matrix
  dac <- ident[,1] # Derive the dac
  w <- ident[,2] # Derive the wmo
  i_nc <- unlist(strsplit(ident[,4],"_"))[2] # Derive the cycle number in the .nc
  #ident[,4] <- i_nc
  n <- as.numeric(strsplit(i_nc,".nc")) # Get the cycle number
  
  # Get the Sfile NetCDF file pathway associated
  file <- list.files(paste(path_data_argo,dac,w,"profiles",sep="/"),
                     pattern = paste("S.*_", formatC(as.integer(n),digits = 2, flag="0"),".nc$",sep=""),full.names = T)
  
  # If file exists
  if(length(file)>0){
    # Get all files for this wmo to know if core variables are in R or D-mode
    # If the DXXXXXX_XXX.nc core file exists then the core parameters are in D-Mode so we load P/T/S data as adjusted data
    # Else the core parameters are in Real time mode
    FILES_IN <- list.files(paste(path_data_argo,dac,w,"profiles",sep="/"),
                           pattern="*.nc") # All files for this WMO
    # Get the mode for physical variables R-mode or D-mode
    # Initialize the Hydro mode in real time
    HYDRO <- "R"
    # If DXXXXXX_XXX.nc core file exists --> hydro mode is D-mode
    if(paste("D",w,"_",paste(formatC(as.integer(n),digits = 2, flag="0"),".nc", sep=""), sep="") %in% FILES_IN){HYDRO <- "D"}
    
    nc <- nc_open(file) # Open the Synthetic netcdf file
    
    # Verify that the Chl is in the file and if not --> continue to next iteration
    STATION_PARAMETERS <- ncvar_get(nc,"STATION_PARAMETERS")
    CHLA_STRING=str_pad("CHLA",64,"right")
    if(!CHLA_STRING %in% STATION_PARAMETERS){
      next
      nc_close(nc)
    }
    
    lon <- ncvar_get(nc,"LONGITUDE") # Get longitude
    lon <- ncvar_get(nc,"LATITUDE") # Get latitude
    # Load core data associated to the HYDRO mode (R or D)
    # If core variables are in D-mode then load the VAR_ADJUSTED parameter
    if(HYDRO=="R"){
      pres <- ncvar_get(nc,"PRES")
      temp <- ncvar_get(nc,"TEMP")
      sal <- ncvar_get(nc,"PSAL")
      pres_qc <- unlist(strsplit(ncvar_get(nc,"PRES_QC"),""))
      temp_qc <- unlist(strsplit(ncvar_get(nc,"TEMP_QC"),""))
      sal_qc <- unlist(strsplit(ncvar_get(nc,"PSAL_QC"),""))
    }
    if(HYDRO=="D"){
      pres <- ncvar_get(nc,"PRES_ADJUSTED")
      temp <- ncvar_get(nc,"TEMP_ADJUSTED")
      sal <- ncvar_get(nc,"PSAL_ADJUSTED")
      pres_qc <- unlist(strsplit(ncvar_get(nc,"PRES_ADJUSTED_QC"),""))
      temp_qc <- unlist(strsplit(ncvar_get(nc,"TEMP_ADJUSTED_QC"),""))
      sal_qc <- unlist(strsplit(ncvar_get(nc,"PSAL_ADJUSTED_QC"),""))
    }
    
    # Remove all bad QC data
    pres[pres_qc %in% c(3,4)] <- NA
    temp[temp_qc %in% c(3,4)] <- NA
    sal[sal_qc %in% c(3,4)] <- NA

#    try({
      
      chl <- ncvar_get(nc,"CHLA_ADJUSTED") # Get chla parameter
      chl_qc <- unlist(strsplit(ncvar_get(nc,"CHLA_ADJUSTED_QC"),""))
      
      nc_close(nc) # Close the NetCDF file
      
      # Remove bad data of CHL and pres
      chl <- chl[!is.na(pres)]
      pres <- pres[!is.na(pres)]
      chl[chl_qc %in% c(3,4)] <- NA
      pres <- pres[!is.na(chl)]
      chl <- chl[!is.na(chl)]
      
      # Condition on chl profile (need profiles...)
      if(length(unique(chl))>5 & !all(is.na(chl)) & !all(is.na(pres)) 
         & length(pres[pres<500])>10 & length(pres[pres>500])>10){
        
        #res <- mean(diff(pres[pres<250]))
        
        # Compute the Chl with new dark correction
        chl_new <- Darkoz(dep_chl = pres, chl = chl)
        
        
        # Plot the profile if no infinite values for chl_new 
        # it happens for example when there is almost no variation in the Chl profile
        if(length(chl_new[!is.infinite(chl_new)])>0){
          plot(chl, -pres, col="black", 
               xlim = range(chl, chl_new, na.rm=T),
               xlab="Chla", ylab = "Depth", 
               main = paste(w, n, sep=" - ")) # Black profile is CHLA_ADJUSTED with normal dark correction
          points(chl_new,-pres, col="red", cex=.5) # Red profiles is CHLA_ADJUSTED with new dark correction
        }
       }
#    })
    print(paste(c, length(chemin_files_to_correct), sep="/"))
  }
}
graphics.off()
