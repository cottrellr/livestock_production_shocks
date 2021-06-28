#Shock detection function based on GAM

#Rich Cottrell 18 June 2021


library(here)
library(tidyverse) 


##Function for shock detection


shock_detector_gam <- function(target_time_series, target_cooks_distance){
  this_time_series <- target_time_series
  this_model <- mgcv::gam(this_time_series~s(seq(1:length(this_time_series))))
  these_residuals <- residuals(this_model)
  these_residuals_t <- these_residuals[2:length(these_residuals)]
  these_residuals_t_1 <- these_residuals[1:length(these_residuals)-1]
  this_regression <- lm(these_residuals_t ~ these_residuals_t_1)
  this_cooks_distance <- cooks.distance(this_regression)
  shock_presence <- c(0, if_else(this_cooks_distance> target_cooks_distance, true = 1, false = 0))
  return(shock_presence)
}

shock_detector_loess <- function(target_time_series, target_cooks_distance, target_span){
  this_time_series <- target_time_series
  this_model <- loess(this_time_series~seq(1:length(this_time_series)), span=target_span)
  these_residuals <- residuals(this_model)
  these_residuals_t <- these_residuals[2:length(these_residuals)]
  these_residuals_t_1 <- these_residuals[1:length(these_residuals)-1]
  this_regression <- lm(these_residuals_t ~ these_residuals_t_1)
  this_cooks_distance <- cooks.distance(this_regression)
  shock_presence <- c(0, if_else(this_cooks_distance> target_cooks_distance, true = 1, false = 0))
  return(shock_presence)
}

#Jessica's function
shockid <- function(dat, thresh){
  # dat is your time series data and threshold is the threshold you want for Cook's D (defaulted to 0.35)
  outt <- array(dim=c(length(dat), 3))
  x <- 1:length(dat)
  ll <- lowess(x, dat) # Fits lowess curve (can specify other options for how the curve is estimated and can change the span)
  rr <- as.numeric(dat[order(x)]-ll$y) #residuals off lowess
  rrp1 <- rr[2:length(rr)] # Residuals at time t
  rrm1 <- rr[1:(length(rr)-1)] # Residuals at time t-1
  
  ll2 <- lm(rrp1~rrm1) # Linear fit of the residuals
  cd <- cooks.distance(ll2) # Calculate the Cook's D
  
  outt[2:length(rr),1] <- as.numeric(cd) # Output the Cook's D
  outt[,2] <- rr # Output the residuals
  outt[2:length(rr),3] <- ifelse(as.numeric(cd) >= thresh,1,0) # Logical of whether point is a shock
  
  return(outt)
}




