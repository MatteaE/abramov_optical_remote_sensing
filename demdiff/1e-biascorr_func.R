# This file contains the definitions for the functions used to examine (and optionally correct)
# various dh biases (across/along track, curvature, altitude, slope, aspect).


#### Function to compute along-track angle of a satellite -----------------------------------------
# This function computes the ground track
# angle of a stereo satellite
# (SPOT 5 or Terra) depending on latitude.
# See Eq. 5.32, Eq. 7.11 and Table 7.1 of
# Capderou M. (2005). Satellites - Orbits and Missions.
# It also handles Pléiades (angle = 0.0, North-South undulations).
# NOTE: orbital_frequency can also be computed by dividing the 1440 minutes in a day
# by the minutes of the orbital period (e.g. for SPOT 5: 1440 / 101.46 = 14 + 5/26).
func_along_track_angle_degrees <- function(latitude_degrees,
                                           scene_type) {
  
  if (scene_type == "SPOT") {
    
    orbital_inclination <- 98.72 * pi/180
    orbital_frequency <- 14 + 5/26
    
  } else if (scene_type == "ASTER") {
    
    orbital_inclination <- 98.21 * pi/180
    orbital_frequency <- 15 - 7/16
    
  } else if (scene_type == "Pléiades") {
    
    return(0.0)
    
  } else {
    
    stop("Scene type for along-track angle calculation must be one of SPOT, ASTER, Pléiades. Received: ", scene_type)
    
  }
  
  i <- orbital_inclination
  k <- orbital_frequency
  phi <- latitude_degrees * pi/180
  
  return(
    atan2(
      cos(i) - (1/k)*cos(phi)*cos(phi),
      sqrt(cos(phi)*cos(phi) - cos(i)*cos(i))
    ) *
      (180 / pi) *
      (-1) # This because we want a positive clockwise angle.
  )
  
}



#### Function to plot a dh raster into bins of an attribute provided in another raster ------------
# This function produces a plot of the (stable-terrain)
# dh provided in dh_ras, classified into bins of the attrib_ras.
# The attribute can be e.g. altitude, curvature, across-track coordinate,
# along-track coordinate.
# class_lim_lower, class_lim_upper, class_step control the binning (and plotting).
func_classify_by_attribute <- function(dh_ras,
                                       attrib_ras,
                                       class_lim_lower,
                                       class_lim_upper,
                                       class_step,
                                       label_x) {
  
  
  dh_v <- values(dh_ras)[,1]
  
  
  # Before the classification, we apply a 4-NMAD filter to dh.
  # This is also what we do when we fit a polynomial or spline to
  # the dh, thus we want to plot the bins in a way that matches the fitted curve.
  # cat(". NMAD filter on dh...\n")
  dh_v[abs(dh_v) > 4 * mad(dh_v, na.rm = TRUE)] <- NA
  
  # cat(". Clamping attribute raster...\n")
  attrib_ras_clamped <- clamp(attrib_ras,
                              lower = class_lim_lower,
                              upper = class_lim_upper,
                              values = FALSE)
  
  
  class_breaks <- seq(class_lim_lower, class_lim_upper, class_step)
  class_midpoints <- (class_breaks[-length(class_breaks)] + class_breaks[-1])/2
  class_rcl_mat <- cbind(class_breaks[-length(class_breaks)],
                         class_breaks[-1],
                         1:length(class_midpoints)) # Classify attribute into integer classes, for safe handling.
  
  # Give some margin at the lowest and highest classes,
  # to avoid floating-point errors in the reclassification
  # (observed while classifying altitude!)
  class_rcl_mat[1,1] <- class_rcl_mat[1,1] - 0.1
  class_rcl_mat[nrow(class_rcl_mat),2] <- class_rcl_mat[nrow(class_rcl_mat),2] + 0.1
  
  # cat(". Classifying attribute raster...\n")
  attrib_class <- classify(attrib_ras_clamped,
                           rcl = class_rcl_mat)
  attrib_class_v <- values(attrib_class)[,1]
  
  
  # Now we need to find the indices of the grid cells belonging to each attribute class.
  # This is very slow if we use a dumb loop with which(... == ...), thus we do it better:
  # we sort the cells by increasing attribute class (first all cells with attribute = 1,
  # then all with = 2, and so on); we keep track of the indices with index.return = TRUE
  # and then we don't need to use which(), we can simply rely on the number of cells within
  # each class (which can be easily computed with rle() since we have sorted by increasing
  # attribute value). It is a real speedup!
  # cat(". Sorting by increasing attribute...\n")
  attrib_sort <- sort(attrib_class_v,
                      index.return = TRUE,
                      na.last = TRUE)
  
  # cat(". RLE...\n")
  attrib_sort_rle <- rle(attrib_sort$x)
  
  # Remove info about NAs, RLE puts each NA in a separate class and we don't care about any of them.
  ids_na <- which(is.na(attrib_sort_rle$values))
  if (length(ids_na) > 0) {
    attrib_sort_rle$lengths <- attrib_sort_rle$lengths[-ids_na]
    attrib_sort_rle$values  <- attrib_sort_rle$values[-ids_na]
  }
  
  # Now deal with the case where some classes are empty (e.g.
  # if the curvature range we consider exceeds the actual
  # curvature values observed in the grid).
  df_out <- data.frame(attrib_class_id = 1:nrow(class_rcl_mat),
                       attrib = class_midpoints,
                       dh = NA_real_,
                       dh_sd = NA_real_,
                       cell_n = NA_integer_)
  
  # These are the indices of all attribute classes within the RLE result.
  # They are NA for each attribute class which is not represented in the input grid.
  df_match_ids <- match(df_out$attrib_class_id, attrib_sort_rle$values) 
  
  # The three below refer ONLY to the classes which are actually represented in the input grid.
  # Thus, later in the loop we select them based on their df_match_ids not being NA.
  attrib_class_sizes <- attrib_sort_rle$lengths
  attrib_class_start_ids <- 1 + c(0, cumsum(attrib_class_sizes)[-nrow(class_rcl_mat)])
  attrib_class_end_ids <- cumsum(attrib_class_sizes)
  if (length(attrib_class_start_ids) > length(attrib_class_end_ids)) { # If not all classes are represented, attrib_class_start_ids has one extra element. We discard it here for tidiness (would not be used anyway).
    attrib_class_start_ids <- attrib_class_start_ids[1:length(attrib_class_end_ids)]
  }
  
  
  # cat(". Computing dh by attribute class...\n")
  for (attrib_class_id in 1:nrow(class_rcl_mat)) {
    # cat("\r", attrib_class_id, "/", nrow(class_rcl_mat))
    
    if (!is.na(df_match_ids[attrib_class_id])) {
      
      ids_cur <- attrib_sort$ix[(attrib_class_start_ids[df_match_ids[attrib_class_id]]):(attrib_class_end_ids[df_match_ids[attrib_class_id]])]
      
      df_out$dh[attrib_class_id] <- mean(dh_v[ids_cur], na.rm = T)
      df_out$dh_sd[attrib_class_id] <- sd(dh_v[ids_cur], na.rm = T)
      df_out$cell_n[attrib_class_id] <- length(which(!is.na(dh_v[ids_cur]))) 
      
    }
  }
  df_out$cell_n[is.na(df_out$cell_n)] <- 0 # NA happens when not only there are no dh values for the current attribute class, but the attribute class itself never appears in the reference attribute grid.
  # cat("\n")
  
  # Plot of the classification.
  pl_out <- ggplot(df_out) +
    geom_tile(aes(x = attrib,
                  fill = cell_n,
                  y = 0),
              width = class_step,
              height = diff(range(c(df_out$dh+df_out$dh_sd, df_out$dh-df_out$dh_sd), na.rm = T))/50) + # Tiles filled with color of cell_n get an automated thickness dependent on the data range.
    geom_hline(yintercept = 0) + #, linewidth = 0.15) +
    geom_errorbar(aes(x = attrib, ymin = dh - dh_sd/2, ymax = dh+dh_sd/2)) + #, linewidth = 0.05) +
    geom_point(aes(x = attrib, y = dh)) + #, size = 0.35, stroke = 0.15) +
    scale_fill_fermenter(palette = "Spectral",
                         limits = c(0,max(df_out$cell_n)),
                         name = "Samples",
                         breaks = 2^pretty(c(0,log2(max(df_out$cell_n))), n = 8)) +
    xlab(label_x) +
    ylab("dh [m]") +
    theme_bw()
  
  return(list(pl_out = pl_out,
              df_out = df_out))
}



#### Function to fit a polynomial or spline correction of stable-terrain dh along an attribute ----
# It returns a list with:
# . numeric vector (same length as ncell(dh_ras)) with the fitted correction for each cell.
# . plot of the binned classification overlain by the fitted correction curve.
# Arguments:
# . dh_ras:              SpatRaster of stable-terrain dh
# . attrib_ras:          SpatRaster with the value of the attribute. Can have NAs (e.g. if it is curvature)
# . attrib_range_fit:    numeric(2), i.e. c(lower,upper). For the fit, disregard values of the attribute outside this range (e.g. extreme curvatures).
# . attrib_class:        numeric(3), i.e. c(lower,upper,step). For the binned plot, disregard values of the attribute outside the given range, and bin according to the given step.
# . attrib_label:        label to use for the attribute, it is put on the X axis of the fit plot
# . plot_ylim:           numeric(2), y (i.e. dh) limits for the plot
# . fit_type:            either "poly<n>" or "spline", where integer <n> is the degree of the polynomial
func_fit_biascorr <- function(dh_ras,
                              attrib_ras,
                              attrib_range_fit,
                              attrib_class,
                              attrib_label,
                              plot_ylim,
                              fit_type) {
  
  
  # Produce plot of the binned dh classification.
  attrib_lower_class <- attrib_class[1]
  attrib_upper_class <- attrib_class[2]
  attrib_step_class  <- attrib_class[3]
  # cat("Plotting binned classification...\n")
  class_list <- func_classify_by_attribute(dh_ras,
                                           attrib_ras,
                                           attrib_lower_class,
                                           attrib_upper_class,
                                           attrib_step_class,
                                           attrib_label)
  pl_class <- class_list$pl_out
  
  
  # Collect stable-terrain dh samples for the fit.
  # We exclude samples where dh is NA (else spline fit fails), and we do 4-NMAD filtering of dh outliers.
  # We also exclude samples with attribute value outside [attrib_lower_fit, attrib_upper_fit].
  # Then we do the fit.
  attrib_lower_fit <- attrib_range_fit[1]
  attrib_upper_fit <- attrib_range_fit[2]
  
  # cat("Collecting samples...\n")
  df_meas_v1 <- data.frame(attrib = values(attrib_ras)[,1],
                           dh = values(dh_ras)[,1])
  df_meas_v2 <- df_meas_v1[which((!is.na(df_meas_v1$dh)) & (!is.na(df_meas_v1$attrib))),]
  df_meas_v3 <- df_meas_v2[which(abs(df_meas_v2$dh) < 4 * mad(df_meas_v1$dh, na.rm = T)),]
  
  rm(df_meas_v2)
  gc()
  
  df_meas_v4 <- df_meas_v3[which((df_meas_v3$attrib >= attrib_lower_fit) & (df_meas_v3$attrib <= attrib_upper_fit)),]
  
  rm(df_meas_v3)
  gc()
  
  # Do the requested fit on the v4 data (the good samples!),
  # apply it (for plotting) over the [attrib_lower_fit, attrib_upper_fit] range.
  # cat("Fitting curve...\n")
  df_plot_mod <- data.frame(attrib = seq(attrib_lower_fit,
                                         attrib_upper_fit,
                                         length.out = 1e4))
  if (fit_type == "spline") {
    
    mod_fit <- smooth.spline(x = df_meas_v4$attrib,
                             y = df_meas_v4$dh,
                             spar = 0.5)
    df_plot_mod$dh <- predict(mod_fit,
                              x = df_plot_mod)$y$attr
    
  } else {
    
    poly_deg <- as.integer(substr(fit_type, 5, nchar(fit_type)))
    mod_fit <- lm(dh~poly(attrib, degree = poly_deg, raw = TRUE),
                  data = df_meas_v4)
    df_plot_mod$dh <- as.numeric(predict.lm(mod_fit,
                                            newdata = df_plot_mod))
    
  }
  
  rm(df_meas_v4)
  gc()
  
  # If the bin classification extends outside the [attrib_lower_fit, attrib_upper_fit] range,
  # we don't extrapolate application of the fitted curve there - rather, we keep the constant
  # value of the edge of the fitting range. In that case, we manually add this value to df_plot_mod.
  if (attrib_lower_class < attrib_lower_fit) {
    df_plot_mod <- rbind(data.frame(attrib = attrib_lower_class,
                                    dh = df_plot_mod$dh[1]),
                         df_plot_mod)
  }
  if (attrib_upper_class > attrib_upper_fit) {
    df_plot_mod <- rbind(df_plot_mod,
                         data.frame(attrib = attrib_upper_class,
                                    dh = df_plot_mod$dh[nrow(df_plot_mod)]))
  }
  
  
  # Extend binned plot with the fitted curve.
  # cat("Plotting fitted curve...\n")
  pl_fit <- pl_class +
    geom_line(data = df_plot_mod,
              aes(x = attrib,
                  y = dh),
              color = "red") +
              # linewidth = 0.2) +
    geom_vline(xintercept = c(attrib_lower_fit, attrib_upper_fit), color = "red") +
    xlim(attrib_lower_class, attrib_upper_class) +
    scale_y_continuous(limits = plot_ylim,
                       oob = scales::oob_keep)
  
  
  # Generate data frame to compute full output correction.
  # We clamp it already with pmin() and pmax(), so that corrections are
  # never extrapolated outside the fitting range.
  # cat("Generating predictor df...\n")
  df_predict <-  data.frame(attrib = pmin(attrib_upper_fit, pmax(attrib_lower_fit, df_meas_v1$attrib)))
  ids_attrib_na <- which(is.na(df_predict$attrib))
  
  rm(df_meas_v1)
  gc()
  
  # The spline fit fails if the predictor attrib is NA,
  # so we temporarily set those to 0.0 and then remove
  # the result again from the corr_vec (i.e., correction = 0.0).
  df_predict$attrib[ids_attrib_na] <- 0.0
  
  # cat("Applying fitted curve...\n")
  if (fit_type == "spline") {
    corr_vec <- predict(mod_fit,
                        x = df_predict)$y$attr
  } else {
    corr_vec <- as.numeric(predict.lm(mod_fit,
                                      newdata = df_predict))
  }
  corr_vec[ids_attrib_na] <- 0.0
  
  return(list(corr_vec = corr_vec,
              pl_fit   = pl_fit))
  
}



#### Function to examine and optionally correct all biases ----------------------------------------
# Biases are: across-track, along-track, max curvature, min curvature, altitude, slope, aspect.
# If the correction is applied, then the next biases are examined (and optionally corrected) on the corrected grid.
# . biascorr_flag is a binary string, e.g. "1210000". 2 means "apply spline correction against this bias", 1 means "apply polynomial correction against this bias", 0 means "examine only, modeling with a polynomial correction".
# . maxc_ras_p and minc_ras_p are paths to the grids of max and min curvature. Those grids should match the origin and resolution of the DEMdiff, so that we can simply use crop() and extend().
# . ele_ref_aras is a SpatRaster used for bias vs altitude, we can just use one of the two DEMs of the DEMdiff.
# . slope_ras_p is the path to the grid of slope, used for bias vs slope. We can recompute and write it on the fly, or use a precomputed one (for Pléiades).
func_biascorr_all <- function(dh_ras_cur,
                              mask_unstable,
                              biascorr_flag,
                              along_track_specification,
                              maxc_ras_p,
                              minc_ras_p,
                              ele_ref_ras,
                              slope_ras_p,
                              aspect_ras_p,
                              plot_dir,
                              plot_ylim_multi) {
  
  biascorr_flags_vec <- as.integer(strsplit(biascorr_flag, split = "")[[1]])
  
  # This and dh_ras_cur will be updated for each bias correction which is actually applied (as opposed to only quantified).
  dh_ras_cur_masked <- mask(dh_ras_cur, mask_unstable, inverse = TRUE) 
  
  
  #### . Across/along track bias ------------------------------------------------------------------
  #### . . Compute shared variables ---------------------------------------------------------------
  ras_x <- setValues(dh_ras_cur, rep(1:ncol(dh_ras_cur), nrow(dh_ras_cur)))
  ras_y <- setValues(dh_ras_cur, rep(1:nrow(dh_ras_cur), each = ncol(dh_ras_cur)))
  
  gc()
  
  scene_center <- rbind(sapply(list(ext(dh_ras_cur)[1:2], ext(dh_ras_cur)[3:4]), mean))
  scene_lat <- project(scene_center, from = crs(dh_ras_cur), to = "EPSG:4326")[,2]
  
  # If we were supplied with the along-track angle (character convertible to numeric),
  # use it. Else compute it from the satellite specification.
  along_track_num <- as.numeric(along_track_specification)
  if (!is.na(along_track_num)) {
    along_track_angle <- along_track_num
  } else {
    along_track_angle <- func_along_track_angle_degrees(latitude_degrees = scene_lat,
                                                        along_track_specification)
  }
  

  
  
  
  #### . . Across-track bias ----------------------------------------------------------------------
  # We do a 5th order polynomial fit to correct across-track stable terrain dh.
  # We don't do spline because there aren't supposed to be any across-track undulations,
  # thus we don't want a high-frequency correction which would be affected by the limited
  # availability of stable terrain (-> random noise).
  cat("Across-track bias... ")
  
  # Number of bins to use for plotting.
  across_bins_n <- 500
  
  ras_across_track_coord <- ras_x * cos(along_track_angle*pi/180) + ras_y * sin(along_track_angle*pi/180)
  
  across_track_minmax <- as.numeric(terra::minmax(ras_across_track_coord))
  across_track_minmax_valid <- as.numeric(terra::minmax(ras_across_track_coord * (dh_ras_cur_masked > -Inf)))
  
  
  across_range_fit   <- across_track_minmax_valid
  across_range_class <- across_track_minmax_valid
  across_class       <- c(across_range_class, (across_range_class[2] - across_range_class[1]) / across_bins_n)
  
  # biascorr_flag values for across-track:
  # 0 = model with poly5, don't apply
  # 1 = model with poly5, apply
  # 2 = model with spline, apply
  if (biascorr_flags_vec[1] %in% c(0,1)) {
    across_mod_type <- "poly5"
  } else if (biascorr_flags_vec[1] == 2) {
    across_mod_type <- "spline"
  } else if (biascorr_flags_vec[1] != 0) {
    stop("Bias correction flag for across-track bias must be one of 0, 1 or 2. Value provided: ", biascorr_flags_vec[1])
  }
  
  across_biascorr_list <- func_fit_biascorr(dh_ras_cur_masked,
                                            ras_across_track_coord,
                                            across_range_fit,
                                            across_class,
                                            "Across-track coordinate",
                                            c(-2, 2) * plot_ylim_multi,
                                            across_mod_type)
  ggsave(file.path(plot_dir, "bins_dh_v2_vs_across.jpg"), across_biascorr_list$pl_fit, width = 20, height = 10)
  
  
  # If asked, apply correction to the demdiff.
  # biascorr_flags_vec[1] can be 1 or 2 and this will be TRUE or 0 and
  # this condition is false so that the correction is not applied.
  if (biascorr_flags_vec[1]) {
    dh_ras_cur <- setValues(dh_ras_cur, values(dh_ras_cur) - across_biascorr_list$corr_vec)
    dh_ras_cur_masked <- mask(dh_ras_cur, mask_unstable, inverse = TRUE)
  }
  
  # Free some memory.
  rm(ras_across_track_coord, across_biascorr_list)
  gc()
  
  
  
  
  #### . . Along-track bias -----------------------------------------------------------------------
  cat("Along-track bias... ")
  
  # Number of bins to use for plotting.
  along_bins_n <- 500
  
  ras_along_track_coord <- -ras_x * sin(along_track_angle*pi/180) + ras_y * cos(along_track_angle*pi/180)
  
  along_track_minmax <- as.numeric(terra::minmax(ras_along_track_coord))
  along_track_minmax_valid <- as.numeric(terra::minmax(ras_along_track_coord * (dh_ras_cur_masked > -Inf)))
  
  
  along_range_fit   <- along_track_minmax_valid
  along_range_class <- along_track_minmax_valid
  along_class       <- c(along_range_class, (along_range_class[2] - along_range_class[1]) / along_bins_n)
  
  # biascorr_flag values for along-track:
  # 0 = model with poly5, don't apply
  # 1 = model with poly5, apply
  # 2 = model with spline, apply
  if (biascorr_flags_vec[2] %in% c(0,1)) {
    along_mod_type <- "poly5"
  } else if (biascorr_flags_vec[2] == 2) {
    along_mod_type <- "spline"
  } else if (biascorr_flags_vec[2] != 0) {
    stop("Bias correction flag for along-track bias must be one of 0, 1 or 2. Value provided: ", biascorr_flags_vec[2])
  }
  
  along_biascorr_list <- func_fit_biascorr(dh_ras_cur_masked,
                                           ras_along_track_coord,
                                           along_range_fit,
                                           along_class,
                                           "Along-track coordinate",
                                           c(-2, 2) * plot_ylim_multi,
                                           along_mod_type)
  ggsave(file.path(plot_dir, "bins_dh_v3_vs_along.jpg"), along_biascorr_list$pl_fit, width = 20, height = 10)
  
  
  # If asked, apply correction to the demdiff.
  # biascorr_flags_vec[2] can be 1 or 2 and this will be TRUE or 0 and
  # this condition is false so that the correction is not applied.
  if (biascorr_flags_vec[2]) {
    dh_ras_cur <- setValues(dh_ras_cur, values(dh_ras_cur) - along_biascorr_list$corr_vec)
    dh_ras_cur_masked <- mask(dh_ras_cur, mask_unstable, inverse = TRUE)
  }
  
  # Free some memory.
  rm(ras_along_track_coord, along_biascorr_list)
  gc()
  
  
  
  #### . Curvature bias ---------------------------------------------------------------------------
  curvature_class_bins_n <- 200 # We plot with this many bins, for both max and min curvature.
  
  # NOTE: we check both maximum and minimum curvature,
  # i.e. the smearing of both ridges and canyons.
  
  #### . . Maximum curvature ----------------------------------------------------------------------
  cat("Max curvature bias... ")
  
  ref_maxc_raw <- rast(maxc_ras_p)
  ref_maxc <- crop(extend(ref_maxc_raw, dh_ras_cur_masked), dh_ras_cur_masked)
  
  if (!compareGeom(ref_maxc, dh_ras_cur_masked, stopOnError = FALSE)) {
    stop("Reference grid of max curvature for bias correction does not match the dh grid!")
  }
  
  # A 4 m Pléiades DEM has a wider curvature range than a 30 m SPOT DEM.
  # We adapt the classification accordingly.
  if (xres(ref_maxc) == 4) {
    maxc_range_fit   <- c(-0.04, 0.20)
    maxc_range_class <- c(-0.20, 0.50)
  } else {
    maxc_range_fit   <- c(-0.002, 0.01)
    maxc_range_class <- c(-0.006, 0.025)
  }
  
  maxc_class         <- c(maxc_range_class, (maxc_range_class[2] - maxc_range_class[1]) / curvature_class_bins_n)
  
  maxc_biascorr_list <- func_fit_biascorr(dh_ras_cur_masked,
                                          ref_maxc,
                                          maxc_range_fit,
                                          maxc_class,
                                          "Maximum curvature",
                                          c(-4, 4) * plot_ylim_multi,
                                          "poly5")
  ggsave(file.path(plot_dir, "bins_dh_v4_vs_maxc.jpg"), maxc_biascorr_list$pl_fit, width = 20, height = 10)
  
  
  # If asked, apply correction to the demdiff.
  if (biascorr_flags_vec[3]) {
    dh_ras_cur <- setValues(dh_ras_cur, values(dh_ras_cur) - maxc_biascorr_list$corr_vec)
    dh_ras_cur_masked <- mask(dh_ras_cur, mask_unstable, inverse = TRUE)
  }
  
  # Free some memory.
  rm(ref_maxc, maxc_biascorr_list)
  gc()
  
  
  
  #### . . Minimum curvature ----------------------------------------------------------------------
  cat("Min curvature bias...  ")
  
  ref_minc_raw <- rast(minc_ras_p)
  ref_minc <- crop(extend(ref_minc_raw, dh_ras_cur_masked), dh_ras_cur_masked)
  
  if (!compareGeom(ref_minc, dh_ras_cur_masked, stopOnError = FALSE)) {
    stop("Reference grid of min curvature for bias correction does not match the dh grid!")
  }
  
  # A 4 m Pléiades DEM has a wider curvature range than a 30 m SPOT DEM.
  # We adapt the classification accordingly.
  if (xres(ref_minc) == 4) {
    minc_range_fit   <- c(-0.20, 0.04)
    minc_range_class <- c(-0.50, 0.20)
  } else {
    minc_range_fit   <- c(-0.01, 0.002)
    minc_range_class <- c(-0.025,  0.006)
  }
  minc_class         <- c(minc_range_class, (minc_range_class[2] - minc_range_class[1]) / curvature_class_bins_n)
  
  minc_biascorr_list <- func_fit_biascorr(dh_ras_cur_masked,
                                          ref_minc,
                                          minc_range_fit,
                                          minc_class,
                                          "Minimum curvature",
                                          c(-4, 4) * plot_ylim_multi,
                                          "poly5")
  ggsave(file.path(plot_dir, "bins_dh_v5_vs_minc.jpg"), minc_biascorr_list$pl_fit, width = 20, height = 10)
  
  
  # If asked, apply correction to the demdiff.
  if (biascorr_flags_vec[4]) {
    dh_ras_cur <- setValues(dh_ras_cur, values(dh_ras_cur) - minc_biascorr_list$corr_vec)
    dh_ras_cur_masked <- mask(dh_ras_cur, mask_unstable, inverse = TRUE)
  }
  
  # Free some memory.
  rm(ref_minc, minc_biascorr_list)
  gc()
  
  
  
  #### . Altitude bias ----------------------------------------------------------------------------
  cat("Altitude bias... ")
  
  if (!compareGeom(ele_ref_ras, dh_ras_cur_masked, stopOnError = FALSE)) {
    stop("Reference grid of altitude for bias correction does not match the dh grid!")
  }
  
  ele_class_n   <- 50
  ele_range_fit <- minmax(ele_ref_ras)
  ele_class     <- c(ele_range_fit, diff(ele_range_fit) / ele_class_n)
  
  ele_biascorr_list <- func_fit_biascorr(dh_ras_cur_masked,
                                         ele_ref_ras,
                                         ele_range_fit,
                                         ele_class,
                                         "Altitude [m]",
                                         c(-2, 2) * plot_ylim_multi,
                                         "poly3")
  ggsave(file.path(plot_dir, "bins_dh_v6_vs_ele.jpg"), ele_biascorr_list$pl_fit, width = 20, height = 10)
  
  
  # If asked, apply correction to the demdiff.
  if (biascorr_flags_vec[5]) {
    dh_ras_cur <- setValues(dh_ras_cur, values(dh_ras_cur) - ele_biascorr_list$corr_vec)
    dh_ras_cur_masked <- mask(dh_ras_cur, mask_unstable, inverse = TRUE)
  }
  
  # Free some memory.
  rm(ele_biascorr_list)
  gc()
  
  
  
  #### . Slope bias -------------------------------------------------------------------------------
  cat("Slope bias...  ")
  
  ref_slope_raw <- rast(slope_ras_p)
  ref_slope <- crop(extend(ref_slope_raw, dh_ras_cur_masked), dh_ras_cur_masked)
  
  if (!compareGeom(ref_slope, dh_ras_cur_masked, stopOnError = FALSE)) {
    stop("Reference grid of slope for bias correction does not match the dh grid!")
  }
  
  slope_range_fit   <- c(0, 90)
  slope_class       <- c(0, 90, 2.5)
  
  slope_biascorr_list <- func_fit_biascorr(dh_ras_cur_masked,
                                           ref_slope,
                                           slope_range_fit,
                                           slope_class,
                                           "Slope [°]",
                                           c(-2, 2) * plot_ylim_multi,
                                           "poly5")
  ggsave(file.path(plot_dir, "bins_dh_v8_vs_slope.jpg"), slope_biascorr_list$pl_fit, width = 20, height = 10)
  
  
  # If asked, apply correction to the demdiff.
  if (biascorr_flags_vec[6]) {
    dh_ras_cur <- setValues(dh_ras_cur, values(dh_ras_cur) - slope_biascorr_list$corr_vec)
    dh_ras_cur_masked <- mask(dh_ras_cur, mask_unstable, inverse = TRUE)
  }
  
  
  # Free some memory.
  rm(ref_slope, slope_biascorr_list)
  gc()
  
  
  
  #### . Aspect bias ------------------------------------------------------------------------------
  cat("Aspect bias... ")
  
  ref_aspect_raw <- rast(aspect_ras_p)
  ref_aspect <- crop(extend(ref_aspect_raw, dh_ras_cur_masked), dh_ras_cur_masked)
  
  if (!compareGeom(ref_aspect, dh_ras_cur_masked, stopOnError = FALSE)) {
    stop("Reference grid of aspect for bias correction does not match the dh grid!")
  }
  
  aspect_range_fit   <- c(0, 360)
  aspect_class       <- c(0, 360, 5)
  
  aspect_biascorr_list <- func_fit_biascorr(dh_ras_cur_masked,
                                            ref_aspect,
                                            aspect_range_fit,
                                            aspect_class,
                                            "Aspect [°]",
                                            c(-0.8, 0.8) * plot_ylim_multi,
                                            "poly5")
  ggsave(file.path(plot_dir, "bins_dh_v8_vs_aspect.jpg"), aspect_biascorr_list$pl_fit, width = 20, height = 10)
  
  
  # If asked, apply correction to the demdiff.
  if (biascorr_flags_vec[7]) {
    dh_ras_cur <- setValues(dh_ras_cur, values(dh_ras_cur) - aspect_biascorr_list$corr_vec)
  }
  
  # Free some memory.
  rm(ref_aspect, aspect_biascorr_list)
  gc()
  
  # This will be different from the input dh_ras_cur only if there is at least a "1" in the corr_flag string.
  return(dh_ras_cur)
  
}

