# This function filters NS and EW velocities with an iterative focal median filter:
# it removes pixels marked as outliers, then recomputes the filter and repeats filtering,
# continuing until convergence (i.e. no more change in the grids).
# Thresholds to mark a pixel as outlier: either median ± n·SD, or median ± n·NMAD,
# or IQR (i.e., Q25 - 1.5·IQR, Q75 + 1.5·IQR).

func_focal_filter_iterate <- function(vel_ew_orig,
                                      vel_ns_orig,
                                      filter_focal_windowsize,
                                      filter_sd_coeff,
                                      iter_max_n) {
  
  # Create these here already, to return them in case the filter does nothing (never happens).
  vel_ew_filter <- vel_ew_orig
  vel_ns_filter <- vel_ns_orig
  
  vel_ew_med <- focal(vel_ew_orig,
                      w = filter_focal_windowsize,
                      median,
                      na.rm = T,
                      na.policy = "all",
                      fillValue = NA)
  vel_ew_sd <- focal(vel_ew_orig,
                     w = filter_focal_windowsize,
                     iqr,
                     na.rm = T,
                     na.policy = "all",
                     fillvalue = NA)
  vel_ns_med <- focal(vel_ns_orig,
                      w = filter_focal_windowsize,
                      median,
                      na.rm = T,
                      na.policy = "all",
                      fillValue = NA)
  vel_ns_sd <- focal(vel_ns_orig,
                     w = filter_focal_windowsize,
                     iqr,
                     na.rm = T,
                     na.policy = "all",
                     fillvalue = NA)
  
  
  # We use a dummy code (9999) to find out how many outliers are found by this filter.
  filter_vel_sd_ns_mask_9999 <- subst(abs(vel_ns_orig - vel_ns_med) < filter_sd_coeff*vel_ns_sd, 0, 9999)
  filter_vel_sd_ew_mask_9999 <- subst(abs(vel_ew_orig - vel_ew_med) < filter_sd_coeff*vel_ew_sd, 0, 9999)

  
  
  filtered_out_n <- length(which((values(filter_vel_sd_ns_mask_9999) == 9999) | (values(filter_vel_sd_ew_mask_9999) == 9999)))
  
  iter_count <- 1
  
  # If we don't enter the loop, we have to apply the filter here.
  if (iter_max_n == 1) {
    filter_vel_sd_ns_mask <- subst(filter_vel_sd_ns_mask_9999, 9999, NA)
    filter_vel_sd_ew_mask <- subst(filter_vel_sd_ew_mask_9999, 9999, NA)
    vel_ns_filter <- vel_ns_orig * filter_vel_sd_ns_mask
    vel_ew_filter <- vel_ew_orig * filter_vel_sd_ew_mask
  }
  
  while ((filtered_out_n > 0) && (iter_count < iter_max_n)) {
    
    cat("Focal filtering iteration", iter_count, "found", filtered_out_n, "outliers.\n")
    iter_count <- iter_count + 1
    
    # If we have actually found outliers, we convert the dummy 9999 into something useful for masking (NA).
    filter_vel_sd_ns_mask <- subst(filter_vel_sd_ns_mask_9999, 9999, NA)
    filter_vel_sd_ew_mask <- subst(filter_vel_sd_ew_mask_9999, 9999, NA)
    vel_ns_filter <- vel_ns_filter * filter_vel_sd_ns_mask
    vel_ew_filter <- vel_ew_filter * filter_vel_sd_ew_mask
    
    vel_ew_med <- focal(vel_ew_filter,
                        w = filter_focal_windowsize,
                        median,
                        na.rm = T,
                        na.policy = "all",
                        fillValue = NA)
    vel_ew_sd <- focal(vel_ew_filter,
                       w = filter_focal_windowsize,
                       iqr,
                       na.rm = T,
                       na.policy = "all",
                       fillvalue = NA)
    vel_ns_med <- focal(vel_ns_filter,
                        w = filter_focal_windowsize,
                        median,
                        na.rm = T,
                        na.policy = "all",
                        fillValue = NA)
    vel_ns_sd <- focal(vel_ns_filter,
                       w = filter_focal_windowsize,
                       iqr,
                       na.rm = T,
                       na.policy = "all",
                       fillvalue = NA)
   
    
    
    filter_vel_sd_ns_mask_9999 <- subst(abs(vel_ns_filter - vel_ns_med) < filter_sd_coeff*vel_ns_sd, 0, 9999)
    filter_vel_sd_ew_mask_9999 <- subst(abs(vel_ew_filter - vel_ew_med) < filter_sd_coeff*vel_ew_sd, 0, 9999)
    filtered_out_n <- length(which((values(filter_vel_sd_ns_mask_9999) == 9999) | (values(filter_vel_sd_ew_mask_9999) == 9999)))
    
  }
  
  
  return(list(vel_ew_filter,
              vel_ns_filter))  
  
}
