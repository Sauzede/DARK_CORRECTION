############## Interpolation function
interp_nan <- function (dep_reg,dep_orig,var_orig) {
  
  if (length(var_orig)>10) {
    var_reg = approx(dep_orig,var_orig,dep_reg)$y
    var_reg[which(dep_reg<min(dep_orig,na.rm = TRUE)+1 | 
                    (dep_reg>max(dep_orig,na.rm = TRUE)-1))]<-NA
  }
  
  else {
    var_reg <- (c(0:max(dep_reg)*NA))
  }
  
  return(var_reg)
  
}