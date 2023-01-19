###############################################################################################
###############################################################################################
######  Function for new dark correction in OMZs and STGs (+ Baltic and Black seas) ###########
###############################################################################################
####################  Raphaëlle Sauzède and Marin Cornec - August 2020 ########################
###############################################################################################
###############################################################################################

# Initial code from Marin Cornec (2019)
# Changed with Raphaëlle Sauzède (August 2020)
# NB: exemple with second peak ok: 3901531

# Function that takes as input the depth of chla measurements, the chla measurements and an index
# The default index is put to 0.5
Darkoz<-function(dep_chl,chl,index=0.5) {
  # Initialize the new chl corrected from the dark
  chl_dark<-chl
  
  #Compute the resolution of measurements
  step_res<-NA
  step_res<-mean(abs(diff(dep_chl[which(dep_chl<=500)])),na.rm = T)
  
  if (step_res!=0) {
    
    # Previous version from Marin Cornec (2019) and Wojtasiewicz et al. (2017): dep_threshold --> 500
    # Version 2020: dep_threshold is the max of the depth
    dep_threshold <- max(dep_chl)
    
    # Interp the pressure to pres_res
    dep_reg_res<-NA
    dep_reg_res <-seq(0,dep_threshold,step_res)
    chl_interp_res = interp_nan(dep_reg_res,dep_chl,chl)
    
    # Filter the chla profile and change the window size filtering according to the resolution of measurements 
    chl_filt_res<-NA
    
    if (step_res < 3) {
      chl_filt_res = RunningFilter(5,chl_interp_res, Method="Median")
    }
    if (step_res > 3) {
      chl_filt_res = RunningFilter(1,chl_interp_res, Method="Median")
    }
    
    # Get the depth of max Chl
    Fchla_max<-NA
    Fchla_max<-dep_reg_res[which.max(chl_filt_res)]
    
    # declare empty variables
    min_chla<-NA
    min_dep<-NA
    
    #identify the depth and the value of min Fchla in the profile
    chl_sub1<-chl_filt_res[which(dep_reg_res>Fchla_max)]
    dep_chl_sub1<-dep_reg_res[which(dep_reg_res>Fchla_max)]
    min_dep<-dep_chl_sub1[min(which.min(chl_sub1))]
    min_chla<-min(chl_sub1,na.rm=T)
    # offset at the minimum Fchla value
    chl_dark<-chl_dark-min_chla
    
    # Test if second peak of Fchla below min_dep
    sub_set_chla<-NA
    sub_set_chla<-chl_filt_res[which(dep_reg_res<=dep_threshold & dep_reg_res>=min_dep)]
    sub_set_dep<-NA
    sub_set_dep<-dep_reg_res[which(dep_reg_res<=dep_threshold & dep_reg_res>=min_dep)]
    
    coef_var<-NA
    coef_var<- sd(sub_set_chla,na.rm=T)/mean(sub_set_chla,na.mr=T)
    
    if(is.na(coef_var)==F) {
      if(coef_var > index) {
        sub_max<-NA
        sub_max<-sub_set_dep[which.max(sub_set_chla)]
        chl_sub2<-NA
        chl_sub2<-chl_filt_res[which(dep_reg_res > sub_max  & dep_reg_res <= dep_threshold)]
        dep_chl_sub2<-NA
        dep_chl_sub2<-dep_reg_res[which(dep_reg_res > sub_max & dep_reg_res <= dep_threshold)]
        
        # Tets if second peak > 0.05
        # bacause of low resolution profiles at depth, when we have noise, it can be considered as a second peak
        if(max(sub_set_chla) > 0.05){min_dep<-dep_chl_sub2[min(which.min(chl_sub2))]}
      }
    }
    chl_dark[which(dep_chl>=min_dep)]<-0
  }
  return(chl_dark)
}