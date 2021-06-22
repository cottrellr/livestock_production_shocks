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
