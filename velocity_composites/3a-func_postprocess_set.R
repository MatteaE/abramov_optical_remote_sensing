# This function takes a folder of correlations (i.e. displacements)
# and returns SpatRasters of median EW and NS annual velocities
# (also performing filtering of the individual displacements)
# as well as the count of grids contributing to each final median pixel.

func_postprocess_set <- function(path_correlations,
                                 filename_format_noext,
                                 dt_sel,
                                 blacklist_dates_str) {
  
  
  #### . Load Cosi-CORR results, compute velocities -------------------------------------------------
  # We load grids found in the output directory, based on their time separation.
  # This allows us to limit the analysis and aggregation only to some pairs,
  # without having to recompute the whole thing.
  
  # List correlations.
  lf_raw <- list.files(path_correlations, pattern = "\\.tif$")
  
  # Discard correlations with blacklisted dates.
  regexp_blacklist <- paste0("(", paste(blacklist_dates_str, collapse = ")|("), ")")
  ids_blacklist <- grep(regexp_blacklist, lf_raw)
  
  cat(length(ids_blacklist), "correlations are discarded using the blacklist.\n")
  
  if (length(ids_blacklist) == length(lf_raw)) {
    stop("blacklist removed all correlations! No mosaic possible!")
  }
  lf <- lf_raw
  if (length(ids_blacklist) > 0) {
    lf <- lf_raw[-ids_blacklist]
  }
  
  
  # Compute number of characters in the original filenames (except for the band specification and file extension).
  # This is a leftover from when we had inconsistent date formats (with or without dashes).
  filename_nchar_noext <- nchar(format(as.Date("2000/01/01"), format = filename_format_noext))
  
  date1_all <- as.Date(substr(lf, 1, filename_nchar_noext), format = filename_format_noext)
  date2_all <- as.Date(substr(lf, filename_nchar_noext + 8, filename_nchar_noext + 8 + filename_nchar_noext - 1), format = filename_format_noext)
  sel_ids <- which((date2_all - date1_all) %in% dt_sel)
  
  lf_sel <- lf[sel_ids]
  cat("Selected", length(lf_sel), "out of", length(lf), "displacement maps, based on the specified time intervals\n")
  
  # Load all grids
  cat("Loading displacement maps...\n")
  
  lp <- file.path(path_correlations, lf_sel)
  displ_all <- rast(lp)
  nl <- nlyr(displ_all)
  displ_maps_n <- nl / 3
  
  # Split bands into EW, NS and SNR.
  displ_ew_rast_all <- displ_all[[seq(1,nl,3)]]
  displ_ns_rast_all <- displ_all[[seq(2,nl,3)]]
  displ_snr_rast_all <- displ_all[[seq(3,nl,3)]]
  
  # Recompute time separation in days only for the selected correlations,
  # taking it from the name of each Cosi-CORR output file.
  fn_all <- basename(sources(displ_ew_rast_all))
  date1_all <- as.Date(substr(fn_all, 1, filename_nchar_noext), format = filename_format_noext)
  date2_all <- as.Date(substr(fn_all, filename_nchar_noext + 8, filename_nchar_noext + 8 + filename_nchar_noext - 1), format = filename_format_noext)
  datediff_all <- as.integer(date2_all - date1_all)
  
  # Convert displacement into annual velocity.
  vel_ew_rast_all <- displ_ew_rast_all * 365.2425 / datediff_all
  vel_ns_rast_all <- displ_ns_rast_all * 365.2425 / datediff_all
  
  
  
  
  #### . Simple filtering ---------------------------------------------------------------------------
  # We compute filtering masks for all grid layers simultaneously,
  # masks have 1 or NA. Then we multiply layer-wise.
  
  cat("Filtering...\n")
  
  # SNR filter
  filter_snr_mask <- subst(displ_snr_rast_all > filter_snr_thresh, 0, NA)
  vel_ew_rast_all_f1 <- vel_ew_rast_all * filter_snr_mask
  vel_ns_rast_all_f1 <- vel_ns_rast_all * filter_snr_mask
  
  
  # Subtract median EW and NS velocities over stable terrain.
  # We do this before magnitude masking, since misregistration
  # of 10 m over 10 days already amounts to 365 m/yr bias and needs
  # to be preserved (not filtered out) in order to subtract it.
  # We use (col)Median(s) and not Mean, since we have not yet filtered for magnitude.
  stable_terrain <- rast(path_stable_terrain)
  stable_terrain_ids <- which(values(stable_terrain) == 1)
  
  vel_ew_mean_stable <- as.numeric(colMedians(as.matrix(vel_ew_rast_all_f1[stable_terrain_ids]), na.rm = T))
  vel_ns_mean_stable <- as.numeric(colMedians(as.matrix(vel_ns_rast_all_f1[stable_terrain_ids]), na.rm = T))
  
  vel_ew_mean_stable[is.na(vel_ew_mean_stable)] <- 0.0
  vel_ns_mean_stable[is.na(vel_ns_mean_stable)] <- 0.0
  
  vel_ew_rast_all_f1 <- vel_ew_rast_all_f1 - vel_ew_mean_stable
  vel_ns_rast_all_f1 <- vel_ns_rast_all_f1 - vel_ns_mean_stable
  
  
  # Velocity magnitude filter.
  vel_mag_f1 <- sqrt(vel_ew_rast_all_f1^2 + vel_ns_rast_all_f1^2)
  filter_vel_mag_mask <- subst(vel_mag_f1 < filter_vel_mag_thresh, 0, NA)
  vel_ew_rast_all_f2 <- vel_ew_rast_all_f1 * filter_vel_mag_mask
  vel_ns_rast_all_f2 <- vel_ns_rast_all_f1 * filter_vel_mag_mask
  vel_mag_rast_all_f2 <- sqrt(vel_ew_rast_all_f2^2 + vel_ns_rast_all_f2^2)
  
  
  # Compute median velocity. Also the unfiltered one, for comparison.
  vel_ew_med_f2 <- terra::median(vel_ew_rast_all_f2, na.rm = T)
  vel_ns_med_f2 <- terra::median(vel_ns_rast_all_f2, na.rm = T)
  
  # Compute velocity magnitude.
  # vel_mag_med <- sqrt(vel_ew_med^2 + vel_ns_med^2)
  vel_mag_med_f2 <- sqrt(vel_ew_med_f2^2 + vel_ns_med_f2^2)
  
  # How many image pairs contributed to the velocity estimation at each pixel?
  # We use the EW set (same result with the NS set).
  # We also compute this relative to the maximum (i.e. to the total number of displacement maps).
  vel_count <- sum(!is.na(vel_ew_rast_all_f2))
  vel_frac <- vel_count / displ_maps_n
  
  
  # #### . Prepare data for velocity plots ------------------------------------------------------------
  # Remove NA collar using terra::trim().
  # vel_mag_med_trim    <- terra::trim(vel_mag_med)
  vel_mag_med_f2_trim <- terra::trim(vel_mag_med_f2)
  vel_ew_med_f2_trim  <- terra::trim(vel_ew_med_f2)
  vel_ns_med_f2_trim  <- terra::trim(vel_ns_med_f2)
  vel_count_trim      <- terra::trim(vel_count, value = 0)  # count is 0 on the collar, we remove it to match the extent of the others.
  vel_frac_trim       <- terra::trim(vel_frac, value = 0.0) # frac is 0.0 on the collar, we remove it.
  
  
  
  #### . Plot all -----------------------------------------------------------------------------------
  cat("Plotting...\n")
  
  pl_vel <- func_plot_velocity(vel_ew_med_f2_trim,
                               vel_ns_med_f2_trim,
                               vel_mag_med_f2_trim)
  
  
  
  return(list(
    vel_ns    = vel_ns_med_f2_trim,
    vel_ew    = vel_ew_med_f2_trim,
    vel_mag   = vel_mag_med_f2_trim,
    vel_count = vel_count_trim,
    vel_frac  = vel_frac_trim,
    pl_vel    = pl_vel
  ))
  
}
