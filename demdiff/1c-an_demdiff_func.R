# This file contains some function definitions used by 1-analyze_demdiff.R.



#### Function to pre-process a DEM (crop, resample, reproject) according to the given flag --------
# 1 = crop/extend to extent, 2 = resample to blueprint, 3 = reproject to blueprint.
func_preprocess <- function(orig_ras_p,
                            blueprint_ras,
                            preproc_flag) {
  orig_ras <- rast(orig_ras_p)
  if (preproc_flag == 1) {
    preproc_ras <- crop(orig_ras, ext(blueprint_ras), extend = TRUE)
  } else if (preproc_flag == 2) {
    preproc_ras <- resample(orig_ras, blueprint_ras, method = "bilinear")
  } else {
    preproc_ras <- project(orig_ras, blueprint_ras, method = "bilinear")
  }
}




#### Function to compute polygon mask of steep pixels ---------------------------------------------
# This function takes a DEM SpatRaster and a slope threshold
# and returns a SpatVector mask with polygons covering all the
# DEM cells whose slope is above the threshold.
func_slope_mask <- function(dem, slope_thresh) {
  
  slope_cur <- terrain(dem,
                       v = "slope",
                       unit = "degrees")
  slope_mask <- slope_cur > slope_thresh
  slope_mask2 <- subst(slope_mask, 0, NA)
  slope_maskv <- as.polygons(slope_mask2)
  return(slope_maskv)
  
}



#### Function to compute polygon mask of 4-NMAD-outliers ------------------------------------------
# This function takes a DEMdiff SpatRaster and returns
# a SpatVector mask with polygons covering all the
# DEMdiff cells whose absolute value is above mad_coeff * NMAD.
# A common choice is mad_coeff = 4.
func_mad_mask <- function(demdiff,
                          mad_coeff) {
  
  mad_cur <- mad(demdiff, na.rm = T) # This is actually the NMAD thanks to default argument "constant = 1.4826"
  dd_mask <- abs(demdiff) > mad_coeff * mad_cur
  dd_mask2 <- subst(dd_mask, 0, NA)
  dd_maskv <- as.polygons(dd_mask2)
  return(dd_maskv)
  
} 



#### Function to compute raster mask of simple threshold outliers ---------------------------------
# This function takes a DEMdiff SpatRaster and returns
# a SpatVector mask with polygons covering all the
# DEMdiff cells whose absolute value is above a threshold.
# It is used to filter out clouds in the ASTER DEMs.
func_dh_mask <- function(demdiff,
                         thresh) {
  
  dd_mask <- abs(demdiff) < thresh
  dd_mask2 <- subst(dd_mask, 0, NA)
  
  return(dd_mask2)
  
}



#### Function to plot a DEM of difference with red/blue colors ------------------------------------
func_plot_demdiff <- function(demdiff_ras,
                              outlines_vec = NULL,
                              colorscale_lim = 10,
                              xlim = ext(demdiff_ras)[1:2],
                              ylim = ext(demdiff_ras)[3:4],
                              theme_base_size = 11) {
  df_plot <- as.data.frame(demdiff_ras, xy=T)
  names(df_plot) <- c("x", "y", "dz")
  pl <- ggplot(df_plot) +
    geom_raster(aes(x = x, y = y, fill = dz)) +
    scale_fill_distiller(palette = "RdBu",
                         direction = 1,
                         limits = c(-colorscale_lim,colorscale_lim),
                         oob = scales::oob_squish,
                         na.value = "#00000000",
                         name = "dh [m]",
                         guide = guide_colourbar(barheight = unit(theme_base_size * theme_base_size / 30, "lines"))) +
    {if (!is.null(outlines_vec)) geom_spatvector(data = outlines_vec, color = "black", fill = NA)} +
    coord_sf(datum = sf::st_crs(32642),
             xlim = xlim,
             ylim = ylim,
             expand = FALSE) +
    scale_x_continuous(name = bquote("UTM X [10"^3*" m]"),
                       labels = function(x) {as.character(x/1e3)}) +
    scale_y_continuous(name = bquote("UTM Y [10"^3*" m]"),
                       labels = function(x) {as.character(x/1e3)}) +
    theme_bw(base_size = theme_base_size) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "#D0D0D0FF"))
  return(pl)
}



#### Function to produce histograms and Q-Q plots on stable ground --------------------------------
# We plot the distribution histogram of the stable-terrain dh.
# We also overlay the corresponding best-fit Gaussian distribution.
# Then we produce quantile-quantile plots.
# We have two versions of the stable terrain dh:
# (1) Original dh, using the stable terrain "v1", as defined for coregistration (i.e. no MAD outliers, no steep slopes and no buffered glaciers)
# (2) Original dh, using the stable terrain "v2", as defined for Hugonnet uncertainty estimation (i.e., only no buffered glaciers).
# NOTE: we do not save the plots yet: we do it later, after running Hugonnet's script,
# and adding to these the plot of standardized uncertainty (which should match a Gaussian quite well).
func_plot_stable_demdiffs <- function(mask_unstable_for_coreg,
                                      demdiff_orig_ras) {
  
  hist_bin_width <- 1 # Width of histogram bins, in meters.
  
  mask_stable_v1_ras <- rasterize(mask_unstable_for_coreg, demdiff_orig_ras, field = NA_real_, background = 1)
  dh_val_stable_v1   <- as.numeric(na.omit(demdiff_orig_ras[mask_stable_v1_ras][,1]))
  mean_v1 <- mean(dh_val_stable_v1)
  sd_v1   <- sd(dh_val_stable_v1)
  
  mask_stable_v2_ras <- rasterize(mask_unstable_for_dh_err, demdiff_orig_ras, field = NA_real_, background = 1)
  dh_val_stable_v2   <- as.numeric(na.omit(demdiff_orig_ras[mask_stable_v2_ras][,1]))
  mean_v2 <- mean(dh_val_stable_v2)
  sd_v2   <- sd(dh_val_stable_v2)
  
  
  pl_dh_stable_v1_hist <- ggplot(data.frame(dh = dh_val_stable_v1)) +
    geom_histogram(aes(x = dh), binwidth = hist_bin_width, center = 0.5) +
    geom_function(fun = function(x){dnorm(x,mean_v1,sd_v1) * hist_bin_width * length(dh_val_stable_v1)}) +
    geom_vline(xintercept = mean_v1, color = "red") +
    xlab("dh over stable terrain (RGIbuf + MAD + slopes) [m]") +
    ylab("Count") +
    theme_bw()
  pl_dh_stable_v1_qq <- ggplot(data.frame(dh = dh_val_stable_v1), aes(sample=dh)) +
    stat_qq() +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Theoretical normal quantile") +
    ylab("Observed quantile") +
    theme_bw()
  
  pl_dh_stable_v2_hist <- ggplot(data.frame(dh = dh_val_stable_v2)) +
    geom_histogram(aes(x = dh), binwidth = hist_bin_width, center = 0.5) +
    geom_function(fun = function(x){dnorm(x,mean_v2,sd_v2) * hist_bin_width * length(dh_val_stable_v2)}) +
    geom_vline(xintercept = mean_v2, color = "red") +
    xlab("dh over stable terrain (RGIbuf only) [m]") +
    ylab("Count") +
    theme_bw()
  pl_dh_stable_v2_qq <- ggplot(data.frame(dh = dh_val_stable_v2), aes(sample=dh)) +
    stat_qq() +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Theoretical normal quantile") +
    ylab("Observed quantile") +
    theme_bw()
  
  
  return(list(pl_dh_stable_v1_hist,
              pl_dh_stable_v1_qq,
              pl_dh_stable_v2_hist,
              pl_dh_stable_v2_qq))
}



#### Function to plot standardized demdiffs (histogram and Q-Q) -----------------------------------
func_plot_stable_demdiffs_std <- function(demdiff_dh_err_std_ras) {
  
  dh_val_standardized <- as.numeric(na.omit(values(demdiff_dh_err_std_ras)))
  
  mean_std <- mean(dh_val_standardized)
  sd_std   <- sd(dh_val_standardized)
  
  pl_dh_stable_std_hist <- ggplot(data.frame(dh = dh_val_standardized)) +
    geom_histogram(aes(x = dh), binwidth = 0.1, center = 0.5) +
    geom_function(fun = function(x){dnorm(x,mean_std,sd_std) * 0.1 * length(dh_val_standardized)}) +
    geom_vline(xintercept = mean_std, color = "red") +
    xlab("Standardized dh over stable terrain (RGIbuf only) [-]") +
    ylab("Count") +
    theme_bw()
  
  pl_dh_stable_std_qq <- ggplot(data.frame(dh = dh_val_standardized), aes(sample=dh)) +
    stat_qq() + 
    geom_abline(slope = 1, intercept = 0) +
    xlab("Theoretical normal quantile") +
    ylab("Observed quantile") +
    theme_bw()
  
  return(list(pl_dh_stable_std_hist,
              pl_dh_stable_std_qq))
  
}



#### Function to apply focal median outlier detection to the DEMdiff grid -------------------------
# NOTE: we don't use this anymore, to find outliers
# we have opted for the 3-SD hypsometric filter below.
func_filter_focal <- function(demdiff_orig_ras) {
  
  focal_ws <- 5
  
  demdiff_focal_med <- focal(demdiff_orig_ras,
                             w = focal_ws,
                             fun = median,
                             na.rm = TRUE,
                             expand = TRUE)
  demdiff_focal_nmad <- focal(demdiff_orig_ras,
                              w = focal_ws,
                              fun = mad,
                              na.rm = TRUE)
  
  demdiff_anom <- abs(demdiff_orig_ras - demdiff_focal_med)
  
  # Fast way to apply the mask: we generate a grid which
  # is NA where the deviation is too big, 0 where it is ok.
  demdiff_filter_ras <- demdiff_orig_ras + (0 / (demdiff_anom < 4*demdiff_focal_nmad))
  
  return(demdiff_filter_ras)
  
  
}



#### Function to apply a local 3-SD filter to the DEMdiff grid ------------------------------------
func_filter_hypso <- function(todo_df,
                              comparison_id,
                              polys_cur,
                              demdiff_orig_ras,
                              nasadem_ref_ras) {
  
  # Select filtering polygons.
  polys_cur_filter_ids_str <- todo_df$outl_ids_filter[comparison_id]
  
  if (polys_cur_filter_ids_str == "all") {
    polys_cur_filter_ids_str <- paste0(polys_cur$id, collapse = "-")
  }
  
  # We never filter on polygons corresponding to transverse profiles or subregions we have already filtered.
  # By removing them here, we can still use "all" in those cases instead of spelling out all polygons.
  polys_cur_filter_ids     <- setdiff(as.integer(strsplit(polys_cur_filter_ids_str, "-", fixed = TRUE)[[1]]), c(400,401,1000,10000,10001,20000))
  polys_cur_filter_all     <- polys_cur[which(polys_cur$id %in% polys_cur_filter_ids)]
  
  # Copy demdiff grid, we filter this one.
  demdiff_filter_ras <- demdiff_orig_ras
  
  # The filter is a 3-SD filter applied locally on each polygon,
  # within 10 elevation bands (but not larger than 50 m span).
  # NOTE: below there is also some code for the NMAD, we used
  # to have a 4-NMAD filter instead. But most published stuff uses 3 · SD.
  polys_cur_filter_n <- length(polys_cur_filter_all)
  for (poly_filter_id_id in 1:polys_cur_filter_n) {
    
    cat("\n ***** Filtering polygon", poly_filter_id_id, "\n")
    
    poly_cur_filter_cur <- polys_cur_filter_all[poly_filter_id_id,]
    
    ele_cur_ref <- terra::extract(nasadem_ref_ras, poly_cur_filter_cur, ID = FALSE, cells = TRUE)
    ele_cur_span <- range(ele_cur_ref[,1], na.rm = T)
    
    # Compute elevation bands: at least 10 bands, not wider than 50 m.    
    eleband_size <- diff(ele_cur_span) / max(10, diff(ele_cur_span) / 50)
    elebands_lower <- seq(ele_cur_span[1],
                          ele_cur_span[2] - eleband_size,
                          eleband_size)
    elebands_higher <- elebands_lower + eleband_size
    elebands <- cbind(elebands_lower,
                      elebands_higher)
    elebands_n <- nrow(elebands)
    elebands[1,1] <- elebands[1,1] - 1 # Give some margin at the lower/upper limits of the extreme bands.
    elebands[elebands_n,2] <- elebands[elebands_n,2] + 1
    
    
    
    # DEMdiff filter: 3·SD, applied on each selected polygon, within
    # elevation bands (at least 10 bands, not wider than 50 m).
    outliers_poly_tot_n <- 0 # Total number of cells filtered out within polygon.
    for (eleband_id in 1:elebands_n) {
      
      # cat("Eleband", eleband_id, "of", elebands_n, paste0("(", round(elebands[eleband_id, 1]), "-", round(elebands[eleband_id, 2]), ")"), "--")
      
      cells_cur_id_id <- which((ele_cur_ref[,1] > elebands[eleband_id, 1]) & (ele_cur_ref[,1] <= elebands[eleband_id, 2]))
      
      # Proceed only if there are cells in the current elevation bands
      # (there should always be some except in case of very very narrow bands
      # arising from tiny filtering polygons).
      if (length(cells_cur_id_id) > 0) {
        
        # cat("", length(cells_cur_id_id), "cells --")
        
        cells_cur_ids <- ele_cur_ref$cell[cells_cur_id_id]
        
        # Find median and SD, and dh values exceeding the outlier threshold.
        dh_val_cur <- demdiff_orig_ras[cells_cur_ids]
        med_cur    <- median(dh_val_cur[,1], na.rm = T)
        mean_cur   <- mean(dh_val_cur[,1], na.rm = T)
        sd_cur     <- sd(dh_val_cur[,1], na.rm = T)
        nmad_cur   <- mad(dh_val_cur[,1], na.rm = T)
        
        outliers_id_id <- which(abs(dh_val_cur[,1] - med_cur) > 3*sd_cur)
        
        # cat(" med/mean/nmad/sd = ", round(med_cur, 1), "/", round(mean_cur, 1), "/", round(nmad_cur, 1), "/", round(sd_cur, 1), "--",
        # "ok range =", round(med_cur - 3*sd_cur, 2), "to", round(med_cur + 3*sd_cur, 2), "--",
        # length(outliers_id_id), "outlier(s)\n")
        outliers_poly_tot_n <- outliers_poly_tot_n + length(outliers_id_id)
        if (length(outliers_id_id) > 0) {
          outliers_id <- cells_cur_ids[outliers_id_id]
          demdiff_filter_ras[outliers_id] <- NA_real_
        } # End if there are outliers.
      } # End if there are cells in the current elevation band.
    } # End iterate on the elevation bands.
    
    cat("Eleband SD filter: removed", outliers_poly_tot_n, "cells\n")
    
    gc()
  } # End iterate on the filtering polygons.
  
  return(demdiff_filter_ras)
  
}



#### Function to plot zoomed-in maps and histograms of the DEMdiff over a polygon -----------------
# We plot both the orig and the filtered versions.
func_plot_poly <- function(demdiff_orig_ras,
                           dh_poly_val_orig,
                           demdiff_filter_ras,
                           dh_poly_val_filter,
                           df_summary,
                           polys_cur_aggr,
                           poly_id_id,
                           colorscale_lim) {
  
  # Plot original DEMdiff, with 5 % margin on all sides of the polygon.
  pl_demdiff_orig_poly <- func_plot_demdiff(demdiff_orig_ras,
                                            outlines_vec = polys_cur_aggr[poly_id_id,],
                                            xlim = ext(polys_cur_aggr[poly_id_id,])[1:2] + 0.05 * c(-diff(ext(polys_cur_aggr[poly_id_id,])[1:2]), diff(ext(polys_cur_aggr[poly_id_id,])[1:2])),
                                            ylim = ext(polys_cur_aggr[poly_id_id,])[3:4] + 0.05 * c(-diff(ext(polys_cur_aggr[poly_id_id,])[3:4]), diff(ext(polys_cur_aggr[poly_id_id,])[3:4])),
                                            colorscale_lim = colorscale_lim)
  
  # Plot distribution histogram of the original dh.
  pl_demdiff_orig_hist_poly <- ggplot(data.frame(dh = dh_poly_val_orig[,1])) +
    geom_histogram(aes(x = dh)) +
    geom_vline(xintercept = df_summary$mean_dh_orig[poly_id_id], color = "red") +
    xlab("Original dh [m]") +
    ylab("Count") +
    theme_bw()
  
  # Plot filtered DEMdiff, with 5 % margin on all sides of the polygon.
  pl_demdiff_filter_poly <- func_plot_demdiff(demdiff_filter_ras,
                                              outlines_vec = polys_cur_aggr[poly_id_id,],
                                              xlim = ext(polys_cur_aggr[poly_id_id,])[1:2] + 0.05 * c(-diff(ext(polys_cur_aggr[poly_id_id,])[1:2]), diff(ext(polys_cur_aggr[poly_id_id,])[1:2])),
                                              ylim = ext(polys_cur_aggr[poly_id_id,])[3:4] + 0.05 * c(-diff(ext(polys_cur_aggr[poly_id_id,])[3:4]), diff(ext(polys_cur_aggr[poly_id_id,])[3:4])),
                                              colorscale_lim = colorscale_lim)
  
  # Plot distribution histogram of the filtered dh.
  pl_demdiff_filter_hist_poly <- ggplot(data.frame(dh = dh_poly_val_filter[,1])) +
    geom_histogram(aes(x = dh)) +
    geom_vline(xintercept = df_summary$mean_dh_filter[poly_id_id], color = "red") +
    xlab("Filtered dh [m]") +
    ylab("Count") +
    theme_bw()
  
  
  # Prepare multi-page PDF plot.
  pl_poly_list <- list(pl_demdiff_orig_poly,
                       pl_demdiff_orig_hist_poly,
                       pl_demdiff_filter_poly,
                       pl_demdiff_filter_hist_poly)
  
  return(pl_poly_list)
}



#### Function to extract dh stats over a polygon and put them in the summary df -------------------
func_compute_polygon_summary <- function(dh_poly_val_filter,
                                         demdiff_dh_err_ras,
                                         df_summary,
                                         poly_id_id,
                                         dir_out_cur_p) {
  
  max_dh_id                             <- which.max(dh_poly_val_filter[,1])
  df_summary$max_dh[poly_id_id]         <- round(dh_poly_val_filter[max_dh_id, 1], 2)
  df_summary$max_dh_err[poly_id_id]     <- round(demdiff_dh_err_ras[dh_poly_val_filter$cell[max_dh_id]][,1], 2)
  max_dh_xy                             <- xyFromCell(demdiff_dh_err_ras, dh_poly_val_filter$cell[max_dh_id])
  df_summary$max_dh_x[poly_id_id]       <- round(max_dh_xy[1], 2)
  df_summary$max_dh_y[poly_id_id]       <- round(max_dh_xy[2], 2)
  
  min_dh_id                             <- which.min(dh_poly_val_filter[,1])
  df_summary$min_dh[poly_id_id]         <- round(dh_poly_val_filter[min_dh_id, 1], 2)
  df_summary$min_dh_err[poly_id_id]     <- round(demdiff_dh_err_ras[dh_poly_val_filter$cell[min_dh_id]][,1], 2)
  min_dh_xy                             <- xyFromCell(demdiff_dh_err_ras, dh_poly_val_filter$cell[min_dh_id])
  df_summary$min_dh_x[poly_id_id]       <- round(min_dh_xy[1], 2)
  df_summary$min_dh_y[poly_id_id]       <- round(min_dh_xy[2], 2)
  
  df_summary$mean_dh_orig[poly_id_id]   <- round(mean(dh_poly_val_orig[,1], na.rm = TRUE), 3)
  df_summary$mean_dh_filter[poly_id_id] <- round(mean(dh_poly_val_filter[,1], na.rm = TRUE), 3)
  
  
  # Each CSV file of Hugonnet uncertainties
  # (produced by script 6b-process_hugonnet_uncertainties.py)
  # has only one line (the result using Gaussian+Spherical
  # variogram model and the best estimate for fully correlated
  # variance. We used to have several lines (one per variogram model)
  # but not anymore. Thus we always select the first (and only) line
  # below to fill df_summary.
  uncertainties_summary_df              <- read.csv(file.path(dir_out_cur_p,
                                                              "hugonnet_uncertainties",
                                                              "integrated_error_poly_plots",
                                                              df_summary$poly_id[poly_id_id],
                                                              "dh_error_integrated.csv"),
                                                    header = T)
  df_summary$mean_dh_err[poly_id_id]   <- uncertainties_summary_df$poly_dh_err_m[1]
  df_summary$fullycorr_pct[poly_id_id] <- uncertainties_summary_df$fullycorr_pct[1]
  
  return(df_summary)
  
}



#### Function to compute "gapfilled" mean dh with local hypsometric mean --------------------------
# This returns a *list* with the hypso mean and the plot of the cubic dh(h) fit.
func_compute_local_hypso_mean <- function(nasadem_ref_ras,
                                          poly_cur,
                                          df_dh_h) {
  
  ele_cur_ref  <- terra::extract(nasadem_ref_ras, poly_cur, ID = FALSE, cells = TRUE)
  ele_cur_span <- range(ele_cur_ref[,1], na.rm = T)
  
  #### . . . . Compute elevation bands: at least 10 bands, not wider than 50 m --------------
  bands_min_n    <- 10
  bands_max_size <- 50
  eleband_size <- diff(ele_cur_span) / max(bands_min_n, diff(ele_cur_span) / bands_max_size)
  elebands_lower <- seq(ele_cur_span[1],
                        ele_cur_span[2] - eleband_size,
                        eleband_size)
  elebands_upper <- elebands_lower + eleband_size
  elebands_df <- data.frame(lower = elebands_lower,
                            upper = elebands_upper,
                            cells_n = 0,
                            mean_dh = NA_real_)
  elebands_n <- nrow(elebands_df)
  elebands_df$lower[1] <- elebands_df$lower[1] - 1 # Give 1 m margin at the lower/upper limits of the extreme bands, for floating-point inaccuracies.
  elebands_df$upper[elebands_n] <- elebands_df$upper[elebands_n] + 1
  
  
  #### . . . . Compute and plot 3rd order dh(h) polynomial ----------------------------------
  # We will use it in case of some empty elevation bands.
  dh_h_mod <- rlm(dh~poly(h,3), data = df_dh_h)
  pl_dh_h_cur <- ggplot(df_dh_h) +
    geom_point(aes(x = h, y = dh), size = 0.1, stroke = 0.1) +
    geom_function(fun = function(x) predict(dh_h_mod, data.frame(h = x)), color = "red") +
    xlab("NASADEM altitude [m]") +
    ylab("dh [m]") +
    coord_flip() +
    theme_bw()
  
  #### . . . . Iterate on the elevation bands -----------------------------------------------
  # For each elevation band, find its cells and compute its size and mean dh.
  for (eleband_id in 1:elebands_n) {
    
    # cat("Eleband", eleband_id, "of", elebands_n, paste0("(", round(elebands_df$lower[eleband_id]), "-", round(elebands_df$upper[eleband_id]), ")"), "--")
    
    cells_cur_id_id <- which((df_dh_h$h > elebands_df$lower[eleband_id]) & (df_dh_h$h <= elebands_df$upper[eleband_id]))
    
    # If there are cells in the current elevation band
    # (there should always be some except in case of very
    # narrow bands arising from tiny aggregation polygons),
    # compute mean dh.
    if (length(cells_cur_id_id) > 0) {
      
      # cat("", length(cells_cur_id_id), "cells --")
      
      elebands_df$cells_n[eleband_id] <- length(cells_cur_id_id)
      elebands_df$mean_dh[eleband_id] <- mean(df_dh_h$dh[cells_cur_id_id], na.rm = TRUE)
    } # End if there are cells in the current elevation band.
    
    # At this stage we have the mean dh of the elevation band,
    # except if the band has no dh data (either because the
    # glacier topography entirely skips the current elevation
    # band (extremely rare) or because the current band is all
    # NA in the DEMdiff within the current aggregation polygon.
    # If that is the case, we compute all dh values of the elevation
    # band using a 3rd degree dh(h) polynomial.
    if (is.na(elebands_df$mean_dh[eleband_id])) {
      
      # cat("no valid dh in cell, using cubic dh(h)...")
      
      band_dh_mod_cur <- predict(dh_h_mod,
                                 data.frame(h = df_dh_h$h[cells_cur_id_id]))
      elebands_df$mean_dh[eleband_id] <- mean(band_dh_mod_cur, na.rm = FALSE)
    } # End computation based on 3rd order polynomial, for empty elevation bands.
    
    # cat("\n")
  } # End iterate on the elevation bands.
  
  # Return computed weighted mean change.
  return(list(pl_dh_h_cur = pl_dh_h_cur,
              mean_dh     = weighted.mean(elebands_df$mean_dh,
                                          w = elebands_df$cells_n)))
  
}
