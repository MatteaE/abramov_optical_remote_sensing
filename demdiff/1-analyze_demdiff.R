# This script computes and analyzes DEMdiffs from selected DEM pairs.

# Inputs for the entire processing are set in 6d-set-input.R. Specifically:
# . A character vector (dem_pool) of DEM paths
# . A character vector (outl_pool) of geo-packages with all the polygons that we consider for dh aggregation. One geo-package for type of comparison (see above); each geo-package can have several polygons depending on the area of interest of each comparison. Field id (unique within each geo-package, usually corresponding to the numbering of our glaciers of interest).
# . A data.frame (todo_df) specifying which comparisons shall take place and their specs.


setwd("./")

#### Load packages --------------------------------------------------------------------------------
library(ggplot2)
library(ggpubr)      # For multi-page PDF plot reports (of each aggregation polygon).
library(imager)
library(lidaRtRee)   # For raster2Cimg() and cimg2Raster()
library(MASS)        # For gap-filling with rlm().
library(matrixStats) # For rowMedians()
library(terra)
library(tidyterra)
source("1c-an_demdiff_func.R")
source("1e-biascorr_func.R")



#### Set options and paths ------------------------------------------------------------------------
# Extents
ext_crop_out       <- ext(700005, 730995, 4377015, 4402995) # This is the extent of all our output grids.
ext_crop_proc      <- ext_crop_out + 300                    # This is the extent that we use for processing (slightly larger to accommodate co-registration).

# Plot options
pl_opt_demdiff_width  <- 5
pl_opt_demdiff_height <- 4

# Paths
dir_out_p                    <- "/PATH/TO/OUTPUT/TREE"    # Here we put the full output, with one subfolder per DEM difference.
dir_preproc_name             <- "preproc/"                # Here we put the preprocessed DEMs (coregistered etc.).
dir_demdiffs_name            <- "demdiffs/"               # Here we put the DEM differences.
dir_demdiffs_biascorr_name   <- "biascorr/"               # Here we put the plots of bias correction/analysis.
dir_demdiffs_poly_name       <- "poly/"                   # Here we put the plots of DEMdiffs over individual polygons.
dir_gapfilling_name          <- "gapfilling_proc/"        # Here we put the plots related to gap-filling.
dir_uncertainties_name       <- "hugonnet_uncertainties/" # Here we put the files and plots related to uncertainty estimation.
mask_unstable_rgibuf_p       <- "/PATH/TO/BUFFERED/OUTLINES/VECTOR/FILE" # Static file with 200 m-buffered RGI 7.0 outlines.
outlines_rgi70_all_p         <- "/PATH/TO/OUTLINES/VECTOR/FILE"          # Static file with all RGI 7.0 outlines over the _proc extent, to plot outlines on top of DEMdiffs.
nasadem_ref_p                <- "/PATH/TO/REFERENCE/DEM/GEOTIFF" # (e.g. NASADEM) reference DEM grid, used to compute topographic parameters for slope-based unstable terrain and for uncertainty quantification, and also for the dh(h) cubic fits.

# Slope threshold above which terrain is considered unstable [°].
slope_threshold_unstable     <- 30

# Thresholds for first quick filter of the DEMdiff.
# Different above and below 4000 m (tongue can change fast, firn area not).
firn_limit        <- 4000 # Height of the separation.
dh_tongue_max     <- 200
dh_firn_max       <- 50



#### Define input files and processings to do -----------------------------------------------------
source("1d-set_input.R")


#### Prepare static data: blueprint SpatRaster and useful GPKGs -----------------------------------
nasadem_ref_larger_ras         <- rast(nasadem_ref_p) # We use this to compute slopes for the unstable terrain mask (larger extent, to avoid edge effects / NAs).
nasadem_ref_orig_ras           <- crop(nasadem_ref_larger_ras, ext_crop_proc)
outlines_rgi70_all_vec         <- vect(outlines_rgi70_all_p) # This is always the same, just for plots (and for threshold filtering on glacier only)



#### Loop on the comparisons to be done (the rows of todo_df) -------------------------------------
todo_n <- nrow(todo_df)

for (comparison_id in 1:todo_n) {
  
  cat("\n", comparison_id, "-- ")
  
  #### . Prepare needed paths ---------------------------------------------------------------------
  tba_orig_p <- dem_pool[todo_df$tba_dem_id[comparison_id]]
  ref_orig_p <- dem_pool[todo_df$ref_dem_id[comparison_id]]
  polys_p    <- outl_pool[todo_df$outl_file_id[comparison_id]]
  
  # Create proc and output subdirectories.
  dir_out_cur_p      <- file.path(dir_out_p, sprintf("%02d", comparison_id))
  dir_preproc_cur_p  <- file.path(dir_out_cur_p, dir_preproc_name)
  dir_demdiffs_cur_p <- file.path(dir_out_cur_p, dir_demdiffs_name)
  dir_biascorr_cur_p <- file.path(dir_demdiffs_cur_p, dir_demdiffs_biascorr_name)
  
  dir.create(dir_out_cur_p)
  dir.create(dir_preproc_cur_p)
  dir.create(dir_demdiffs_cur_p)
  dir.create(dir_biascorr_cur_p)
  
  
  #### . Prepare unstable mask used for uncertainty -----------------------------------------------
  mask_unstable_buf_vec      <- vect(unstable_outlines_pool[todo_df$coreg_outl_id[comparison_id]])
  
  mask_unstable_for_dh_err   <- mask_unstable_buf_vec
  mask_unstable_for_dh_err   <- rbind(mask_unstable_for_dh_err, func_slope_mask(nasadem_ref_larger_ras, 45))
  
  mask_unstable_for_dh_err_p <- file.path(dir_preproc_cur_p, "mask_unstable_for_dh_err.gpkg")
  writeVector(mask_unstable_for_dh_err, mask_unstable_for_dh_err_p, overwrite = TRUE, options = "")
  
  
  
  #### . Pre-process DEMs (cropping, resampling, reprojecting) --------------------------------------
  
  preproc_blueprint <- nasadem_ref_orig_ras
  polys_cur         <- vect(polys_p) # Aggregation outlines
  
  # Masks for threshold filtering of on-glacier DEMdiff:
  # mask_firn_v is everything above firn_limit (only within RGI and polygons),
  # mask_tongue_v is the same but below firn_limit.
  mask_firn_v       <- as.polygons(subst(preproc_blueprint >= firn_limit, 0, NA))
  mask_tongue_v     <- as.polygons(subst(preproc_blueprint < firn_limit, 0, NA))
  
  
  cat("Pre-processing DEMs...\n")
  
  # Pre-process TBA DEM.
  tba_preproc_ras <- func_preprocess(tba_orig_p,
                                     preproc_blueprint,
                                     todo_df$tba_dem_preproc[comparison_id])
  tba_preproc_p   <- file.path(dir_preproc_cur_p, "dem_tba_preproc.tif")
  writeRaster(tba_preproc_ras, tba_preproc_p, overwrite = TRUE)
  
  # Pre-process REF DEM.
  ref_preproc_ras <- func_preprocess(ref_orig_p,
                                     preproc_blueprint,
                                     todo_df$ref_dem_preproc[comparison_id])
  ref_preproc_p   <- file.path(dir_preproc_cur_p, "dem_ref_preproc.tif")
  writeRaster(ref_preproc_ras, ref_preproc_p, overwrite = TRUE)
  
  
  
  #### . If needed, co-register TBA to REF ----------------------------------------------------------
  # We generate the unstable terrain mask for coregistration in any case (even if we
  # don't need to do coregistration), just to examine the distribution histogram of dh over
  # stable terrain also for this (stricter) definition of stable terrain.
  # For coregistration, we define unstable ground as the union of:
  # . RGI7.0 outlines, extended with a 200 m buffer
  # . All NASADEM slopes greater than a threshold
  # . All points for which (before coregistration) abs(dh) is > 4·NMAD(dh)
  # The first two were loaded from a static file.
  
  cat("Generating unstable terrain mask...\n")
  demdiff_nocoreg_ras <- tba_preproc_ras - ref_preproc_ras
  mask_unstable_mad_vec      <- func_mad_mask(demdiff_nocoreg_ras, mad_coeff = 4)
  mask_unstable_slope_vec    <- func_slope_mask(nasadem_ref_larger_ras, slope_threshold_unstable)
  mask_unstable_for_coreg    <- rbind(mask_unstable_buf_vec, mask_unstable_mad_vec, mask_unstable_slope_vec)
  
  
  demdiff_nocoreg_masked <- mask(demdiff_nocoreg_ras, mask_unstable_for_coreg, inverse = TRUE)
  stable_frac <- length(which(!is.na(values(demdiff_nocoreg_masked)[,1]))) / length(which(!is.na(values(demdiff_nocoreg_ras)[,1])))
  cat("Available stable terrain:", round(100 * stable_frac, 1), "%")
  
  
  if (todo_df$tba_coreg_flag[comparison_id] > 0) {
    
    #### . . Produce gpkg of unstable ground -----------------------------------------------------------
    mask_unstable_for_coreg_p  <- file.path(dir_preproc_cur_p, "mask_unstable_for_coreg.gpkg")
    writeVector(mask_unstable_for_coreg, mask_unstable_for_coreg_p, overwrite = TRUE, options = "")
    
    
    cat("\nRunning co-registration...\n")
    tba_coreg_out_p <- file.path(dir_preproc_cur_p, "dem_tba_coreg.tif")
    cmd <- "python"
    args <- c("./1a-coreg_local.py",
              ref_preproc_p,
              tba_preproc_p,
              mask_unstable_for_coreg_p,
              tba_coreg_out_p)
    system2(cmd, args)
    
    
    
    # Else: skip coregistration, just use pre-processed TBA as DEM to compute dh.
  } else {
    
    cat("\nSkipping co-registration...\n")
    tba_coreg_out_p <- tba_preproc_p
    
    # End "else, skip coregistration".
  }
  
  # Load coregistered TBA DEM.
  tba_coreg_ras    <- rast(tba_coreg_out_p)
  
  
  
  #### . Compute the DEMdiff ----------------------------------------------------------------------
  # We also plot the DEMdiff histogram over stable terrain
  # (defined as per Hugonnet's script: only RGIbuf).
  demdiff_orig_ras <- tba_coreg_ras - ref_preproc_ras
  
  
  
  #### . First quick filter (DEMdiff threshold) ---------------------------------------------------
  # The threshold is applied separately above and below the firn line.
  mask_firn_r   <- rasterize(mask_firn_v, demdiff_orig_ras, background = 0)
  mask_tongue_r <- rasterize(mask_tongue_v, demdiff_orig_ras, background = 0)

  demdiff_orig_ras[abs(demdiff_orig_ras * mask_firn_r) > dh_firn_max] <- NA_real_
  demdiff_orig_ras[abs(demdiff_orig_ras * mask_tongue_r) > dh_tongue_max] <- NA_real_

  
  
  #### . Store and plot the DEMdiff before bias correction ----------------------------------------
  # We do this even if we don't later apply the bias correction.
  demdiff_prebiascorr_p   <- file.path(dir_demdiffs_cur_p, "demdiff_prebiascorr.tif")
  writeRaster(demdiff_orig_ras, demdiff_prebiascorr_p, overwrite = TRUE)
  
  pl_demdiff_orig <- func_plot_demdiff(demdiff_orig_ras,
                                       outlines_vec = outlines_rgi70_all_vec,
                                       colorscale_lim = 10 * todo_df$plot_scale_mult[comparison_id])
  ggsave(filename = file.path(dir_demdiffs_cur_p, "demdiff_prebiascorr_plot.pdf"),
         plot = pl_demdiff_orig,
         width = pl_opt_demdiff_width,
         height = pl_opt_demdiff_height)
  
  
  
  #### . Estimate and if needed correct the DEMdiff biases ----------------------------------------
  # We quantify and model the following biases:
  # . Across-track (5th order polynomial)
  # . Along-track (spline)
  # . Max and min curvature (5th order polynomial)
  # . Altitude
  # . Slope
  # . Aspect
  # For each, we plot the dh distribution on stable terrain.
  # Afterwards, we will check whether any of those are actually worth correcting.
  # NOTE: we keep the name demdiff_orig_ras, i.e., our original DEMdiff is the bias-corrected one.
  cat("Estimating stable-terrain dh biases... ")
  demdiff_postbiascorr_ras <- func_biascorr_all(demdiff_orig_ras,
                                                mask_unstable_for_coreg,
                                                todo_df$biascorr_flag[comparison_id],
                                                todo_df$biascorr_along_track[comparison_id],
                                                maxc_pool[todo_df$biascorr_topogrids_id[comparison_id]],
                                                minc_pool[todo_df$biascorr_topogrids_id[comparison_id]],
                                                ref_preproc_ras,
                                                slope_pool[todo_df$biascorr_topogrids_id[comparison_id]],
                                                aspect_pool[todo_df$biascorr_topogrids_id[comparison_id]],
                                                dir_biascorr_cur_p,
                                                bias_plots_ylim_multi_pool[todo_df$biascorr_topogrids_id[comparison_id]])
  cat("Done.\n")
  gc()
  
  
  
  #### . If needed, store and plot the bias-corrected DEMdiff -------------------------------------
  # We do this only if we have applied the bias correction.
  if (todo_df$biascorr_flag[comparison_id] != "0000000") {
    
    demdiff_postbiascorr_p   <- file.path(dir_demdiffs_cur_p, "demdiff_postbiascorr.tif")
    writeRaster(demdiff_postbiascorr_ras, demdiff_postbiascorr_p, overwrite = TRUE)
    
    pl_demdiff_postbiascorr <- func_plot_demdiff(demdiff_postbiascorr_ras,
                                                 outlines_vec = outlines_rgi70_all_vec,
                                                 colorscale_lim = 10 * todo_df$plot_scale_mult[comparison_id])
    ggsave(filename = file.path(dir_demdiffs_cur_p, "demdiff_postbiascorr_plot.pdf"),
           plot = pl_demdiff_postbiascorr,
           width = pl_opt_demdiff_width,
           height = pl_opt_demdiff_height)
  }
  
  
  
  # This will stay like this if we don't coregister
  # again, or become the coregistered one if we do.
  demdiff_cur_ras <- demdiff_postbiascorr_ras
  
  
  
  #### . Optional second coregistration after bias correction -------------------------------------
  if (todo_df$tba_coreg_flag[comparison_id] == 2) {
    
    tba_postbiascorr_ras <- ref_preproc_ras + demdiff_postbiascorr_ras
    tba_postbiascorr_p <- file.path(dir_demdiffs_cur_p, "dem_tba_postbiascorr.tif")
    writeRaster(tba_postbiascorr_ras, tba_postbiascorr_p, overwrite = TRUE)
    
    cat("\nRunning second co-registration...\n")
    tba_coreg_out_p <- file.path(dir_preproc_cur_p, "dem_tba_postbiascorr_coreg.tif")
    cmd <- "python"
    args <- c("./1a-coreg_local.py",
              ref_preproc_p,
              tba_postbiascorr_p,
              mask_unstable_for_coreg_p,
              tba_coreg_out_p)
    system2(cmd, args)
    demdiff_postbiascorr_coreg_ras <- rast(tba_coreg_out_p) - ref_preproc_ras
    
    
    
    #### . . Store and plot DEMdiff after second coregistration -----------------------------------
    demdiff_postbiascorr_coreg_p   <- file.path(dir_demdiffs_cur_p, "demdiff_postbiascorr_coreg.tif")
    writeRaster(demdiff_postbiascorr_coreg_ras, demdiff_postbiascorr_coreg_p, overwrite = TRUE)
    
    pl_demdiff_orig <- func_plot_demdiff(demdiff_postbiascorr_coreg_ras,
                                         outlines_vec = outlines_rgi70_all_vec,
                                         colorscale_lim = 10 * todo_df$plot_scale_mult[comparison_id])
    ggsave(filename = file.path(dir_demdiffs_cur_p, "demdiff_postbiascorr_coreg_plot.pdf"),
           plot = pl_demdiff_orig,
           width = pl_opt_demdiff_width,
           height = pl_opt_demdiff_height)
    
    # Update current DEMdiff in case we have coregistered again.
    demdiff_cur_ras <- demdiff_postbiascorr_coreg_ras
    
    rm(demdiff_postbiascorr_coreg_ras, tba_postbiascorr_ras, demdiff_postbiascorr_ras)
    gc()
  }
  
  
  # The DEMdiff for uncertainty calculation is the one after
  # all possible processing (i.e. including 2nd coreg if done).
  demdiff_for_uncertainty_p <- file.path(dir_demdiffs_cur_p, "demdiff_for_uncertainty.tif")
  writeRaster(demdiff_cur_ras, demdiff_for_uncertainty_p, overwrite = TRUE)
  
  
  
  #### . Plot DEMdiff over stable terrain (histogram and Q-Q) -------------------------------------
  pl_stable_terrain_list <- func_plot_stable_demdiffs(mask_unstable_for_coreg,
                                                      demdiff_cur_ras)
  
  
  
  
  #### . Select aggregation outlines --------------------------------------------------------------
  polys_cur_aggr_ids_str  <- todo_df$outl_ids_aggr[comparison_id]
  
  if (polys_cur_aggr_ids_str == "all") {
    polys_cur_aggr_ids_str <- paste0(polys_cur$id, collapse = "-")
  }
  
  polys_cur_aggr_ids     <- as.integer(strsplit(polys_cur_aggr_ids_str, "-", fixed = TRUE)[[1]])
  polys_cur_aggr         <- polys_cur[which(polys_cur$id %in% polys_cur_aggr_ids)]
  
  
  
  #### . Run Hugonnet uncertainty quantification ----------------------------------------------------
  # NOTE: we do this before filtering the DEMdiff grid,
  # anyway we later filter only on the glacier polygons,
  # which are not used to estimate uncertainty (only stable
  # terrain is used). Thus in the end the result is the same.
  cat("Computing uncertainties...\n")
  dir_uncertainties_p <- file.path(dir_out_cur_p, dir_uncertainties_name)
  dir.create(dir_uncertainties_p)
  
  env_mpl <- "MPLBACKEND=agg" # Needed to avoid core dump due to threading issue.
  
  cmd <- "python"
  args <- c("./1b-process_hugonnet_uncertainty.py",
            demdiff_for_uncertainty_p,
            nasadem_ref_p,
            polys_p,
            polys_cur_aggr_ids_str,
            mask_unstable_for_dh_err_p,
            dir_uncertainties_p)
  system(paste(env_mpl, cmd, paste0(args, collapse = " "), sep = " "))
  
  
  # Load raster of dh uncertainties, we will sample it later
  # for the uncertainty of the min/max changes.
  demdiff_dh_err_p <- file.path(dir_uncertainties_p, "dh_error.tif")
  demdiff_dh_err_ras <- rast(demdiff_dh_err_p)
  
  
  #### . . Plot standardized dh over stable terrain (histogram and Q-Q) ---------------------------
  # These should match a Gaussian distribution quite closely.
  # Load raster of standardized dh uncertainties.
  demdiff_dh_err_std_p <- file.path(dir_uncertainties_p, "dh_error_standardized.tif")
  demdiff_dh_err_std_ras <- rast(demdiff_dh_err_std_p)
  
  pl_stable_terrain_list <- c(pl_stable_terrain_list,
                              func_plot_stable_demdiffs_std(demdiff_dh_err_std_ras))
  
  # Export histograms and Q-Q plots of dh over stable terrain (as in Hugonnet22 Fig. S10).
  cat("Producing stable-terrain histograms and Q-Q plots...\n")
  ggexport(list(egg::ggarrange(plots = pl_stable_terrain_list[1:2], ncol = 1, draw = FALSE),
                egg::ggarrange(plots = pl_stable_terrain_list[3:4], ncol = 1, draw = FALSE),
                egg::ggarrange(plots = pl_stable_terrain_list[5:6], ncol = 1, draw = FALSE)),
           filename = file.path(dir_demdiffs_cur_p, "demdiff_stable_hist.pdf"),
           width = 5,
           height = 8)
  
  
  
  
  #### . Filter demdiff over the specified polygons -----------------------------------------------
  cat("Filtering the DEMdiff grid...\n")
  
  # We filter with a 3·SD polygon-wise hypsometric filter,
  # For the Pléiades-Pléiades 4 m comparison, we have to upsample the
  # reference DEM to be able to perform filtering and gapfilling.
  # We do it only once here (not once per each polygon!).
  # We redefine nasadem_ref_ras at each comparison since we may use it at higher
  # resolution for Pléiades and then again the original version for another comparison.
  nasadem_ref_ras <- nasadem_ref_orig_ras
  if (xres(demdiff_cur_ras) != xres(nasadem_ref_orig_ras)) {
    cat("Resampling reference NASADEM to the resolution of the DEMdiff...\n")
    nasadem_ref_ras <- resample(nasadem_ref_orig_ras, demdiff_cur_ras, method = "bilinear")
  }
  demdiff_filter_ras <- func_filter_hypso(todo_df,
                                          comparison_id,
                                          polys_cur,
                                          demdiff_cur_ras,
                                          nasadem_ref_ras)
  
  
  #### . . Save and plot filtered demdiff grid ----------------------------------------------------
  demdiff_filter_p   <- file.path(dir_demdiffs_cur_p, "demdiff_filter.tif")
  writeRaster(demdiff_filter_ras, demdiff_filter_p, overwrite = TRUE)
  
  pl_demdiff_filter <- func_plot_demdiff(demdiff_filter_ras,
                                         outlines_vec = outlines_rgi70_all_vec,
                                         colorscale_lim = 10 * todo_df$plot_scale_mult[comparison_id])
  ggsave(filename = file.path(dir_demdiffs_cur_p, "demdiff_filter_plot.pdf"),
         plot = pl_demdiff_filter,
         width = pl_opt_demdiff_width,
         height = pl_opt_demdiff_height)
  
  
  
  #### . Declare output df for polygon aggregation ------------------------------------------------
  df_summary <- data.frame(
    poly_id             = rep(NA_integer_,
                              nrow(polys_cur_aggr)),        # Id of the aggregation polygon within the GPKG file
    poly_area           = NA_real_,                         # Area of the polygon, computed from the dh cells selected by the polygon
    dh_valid_orig_pct   = NA_real_,                         # Percent of valid dh cells within the polygon (original grid)
    dh_valid_filter_pct = NA_real_,                         # Percent of valid dh cells within the polygon (filtered grid)
    dh_filtered_n       = NA_integer_,                      # Number of filtered cells
    max_dh              = NA_real_,                         # Max dh within each polygon (all max and min are on the original dh grid)
    max_dh_err          = NA_real_,                         # Uncertainty of this max dh, sampled from Hugonnet heteroscedasticity map
    max_dh_x            = NA_real_,                         # X coordinate of the max
    max_dh_y            = NA_real_,                         # Y coordinate of the max
    min_dh              = NA_real_,                         # Min dh within each polygon
    min_dh_err          = NA_real_,                         # Uncertainty of it
    min_dh_x            = NA_real_,                         # X coordinate of the min
    min_dh_y            = NA_real_,                         # Y coordinate of the min
    mean_dh_orig        = NA_real_,                         # Mean dh within each polygon, using not-filtered dh grid
    mean_dh_filter      = NA_real_,                         # Mean dh within each polygon, using filtered dh grid
    mean_dh_hypso       = NA_real_,                         # Mean dh within each polygon, computed with hypsometric local mean as in McNabb2019
    mean_dh_err         = NA_real_,                         # Uncertainty of mean dh within each polygon, computed by Hugonnet script from double-range exponential variogram
    dvol_hypso          = NA_real_,                         # Volume change within each polygon, computed from mean_dh_hypso
    dvol_hypso_err      = NA_real_,                         # Error in the previous one
    fullycorr_pct       = NA_real_                          # Percentage of fully correlated variance needed to best match variogram integration with patches method
  )
  
  
  
  
  #### . Loop over the polygons of aggregation ----------------------------------------------------
  # In this section we do the following:
  # . Compute aggregate values over each polygon
  # . Compute aggregate dh with hypsometric local mean method.
  # . Plot zoomed-in dh grids (orig and filtered) with just the polygon,
  #   and distribution histogram of original and filtered dh grids.
  polys_n <- nrow(polys_cur_aggr)
  polys_size <- expanse(polys_cur_aggr)
  
  dir_demdiffs_poly_p <- file.path(dir_demdiffs_cur_p, "poly")
  dir.create(dir_demdiffs_poly_p)
  
  
  #### . Iterate on the aggregation polygons ------------------------------------------------------
  for (poly_id_id in 1:polys_n) {  # poly_id_id is the row number of the polygon in our loaded GPKG file.
    
    poly_id <- polys_cur_aggr$id[poly_id_id]
    
    cat("Polygon", poly_id, paste0("(", poly_id_id, "/", polys_n, ")"), "\n")
    
    poly_cur <- polys_cur_aggr[poly_id_id,]
    
    #### . . Extract summary statistics of the current polygon --------------------------------------
    dh_poly_val_orig       <- terra::extract(demdiff_cur_ras, poly_cur, ID = FALSE, cells = TRUE)
    dh_poly_val_filter     <- terra::extract(demdiff_filter_ras, poly_cur, ID = FALSE, cells = TRUE)
    
    dh_poly_na_orig_ids    <- which(is.na(dh_poly_val_orig[,1])) # NOTE: these are indices relative to the raster subset covered by the polygon. The absolute cell numbers are in dh_poly_val_orig$cell.
    dh_poly_na_filter_ids  <- which(is.na(dh_poly_val_filter[,1])) # NOTE: these are indices relative to the raster subset covered by the polygon. The absolute cell numbers are in dh_poly_val_orig$cell.
    
    dh_poly_na_orig_n      <- length(dh_poly_na_orig_ids)
    dh_poly_na_filter_n    <- length(dh_poly_na_filter_ids)
    
    df_summary$poly_id[poly_id_id] <- poly_id
    df_summary$poly_area[poly_id_id] <- nrow(dh_poly_val_orig) * xres(demdiff_cur_ras) * yres(demdiff_cur_ras)
    
    df_summary$dh_valid_orig_pct[poly_id_id] <- round(100 * (1 - (dh_poly_na_orig_n / nrow(dh_poly_val_orig))),2)
    df_summary$dh_valid_filter_pct[poly_id_id] <- round(100 * (1 - (dh_poly_na_filter_n / nrow(dh_poly_val_filter))),2)
    df_summary$dh_filtered_n[poly_id_id] <- dh_poly_na_filter_n - dh_poly_na_orig_n
    
    # If polygon has any dh data after filtering:
    # . Find max and min change (over the filtered data), their locations and uncertainties
    # . Find mean dh change (both original and filtered)
    # . Compute mean "gapfilled" dh change with local hypsometric mean
    # . Plot polygon: one PDF with (1) orig dh map, (2) orig dh histogram within polygon, (3) filtered dh map, (4) filtered dh histogram within polygon
    if (df_summary$dh_valid_filter_pct[poly_id_id] > 0) {
      
      # This fills in max/min/mean dh.
      df_summary <- func_compute_polygon_summary(dh_poly_val_filter,
                                                 demdiff_dh_err_ras,
                                                 df_summary,
                                                 poly_id_id,
                                                 dir_out_cur_p)
      # Update mean change error for data gaps, using a conservative factor of five
      # for their uncertainty (Berthier et al., 2014; Dussaillant et al., 2018, Eq. 1)
      df_summary$mean_dh_err[poly_id_id] <- df_summary$mean_dh_err[poly_id_id] * (5 * (100 - df_summary$dh_valid_filter_pct[poly_id_id])/100 + df_summary$dh_valid_filter_pct[poly_id_id]/100)
      
      
      #### . . Compute "gapfilled" mean dh using local hypsometric mean ---------------------------
      # We split the dh in elevation bands, at least 10 and not wider than 50 m.
      # Then we compute the mean change for each band and the final mean change
      # as area-weighted mean of the bands.
      # In case a band has no valid dh data, we fill it with data from a 3rd
      # order dh(h) fit on the whole-polygon data.
      df_dh_h <- data.frame(cell = dh_poly_val_filter$cell,
                            dh   = dh_poly_val_filter[,1],
                            h    = nasadem_ref_ras[dh_poly_val_filter$cell][,1])
      
      #### . . . Proceed only if the polygon has at least 50 % of valid dh values -----------------
      dh_valid_n <- length(which(!is.na(df_dh_h$dh)))
      if (dh_valid_n > 0.5 * nrow(df_dh_h)) {
        
        hypso_mean_list <- func_compute_local_hypso_mean(nasadem_ref_ras,
                                                         poly_cur,
                                                         df_dh_h)
        
        df_summary$mean_dh_hypso[poly_id_id]  <- round(hypso_mean_list$mean_dh, 3)
        df_summary$dvol_hypso[poly_id_id]     <- round(df_summary$poly_area[poly_id_id] * df_summary$mean_dh_hypso[poly_id_id])
        df_summary$dvol_hypso_err[poly_id_id] <- round(df_summary$poly_area[poly_id_id] * df_summary$mean_dh_err[poly_id_id])
        
      } else { # End compute local hypsometric mean only if the polygon has at least 50 % valid dh values.
        
        cat("Warning: polygon has less than 50 % data coverage, I am skipping calculation of hypsometric mean dh.\n")
        
      }
      
      
      #### . . Plot zoomed-in maps and histograms of polygon dh -----------------------------------
      pl_poly_list <- func_plot_poly(demdiff_cur_ras,
                                     dh_poly_val_orig,
                                     demdiff_filter_ras,
                                     dh_poly_val_filter,
                                     df_summary,
                                     polys_cur_aggr,
                                     poly_id_id,
                                     colorscale_lim = 10 * todo_df$plot_scale_mult[comparison_id])
      # Export PDF plots of polygon dh (maps and histograms, original and filtered).
      # We also append the plot of the 3rd order dh(h) polynomial (even in case we didn't use it).
      poly_plots_cur_p <- file.path(dir_demdiffs_poly_p, paste0(poly_id, "_plots.pdf"))
      
      if (dh_valid_n > 0.5 * nrow(df_dh_h)) {
        pl_poly_list <- c(pl_poly_list, hypso_mean_list["pl_dh_h_cur"])
      }
      
      ggexport(pl_poly_list,
               filename = poly_plots_cur_p,
               width = pl_opt_demdiff_width,
               height = pl_opt_demdiff_height)
      
    } # End process polygon if it has dh data.
  } # End of loop over the polygons of interest.
  
  
  
  #### . Write df_summary to file -----------------------------------------------------------------
  # We sort by increasing poly_id, thus
  # Abramov will always be among the first lines.
  ids_sort <- sort(df_summary$poly_id, index.return = TRUE)
  df_summary_out <- df_summary[ids_sort$ix,]
  
  write.csv(df_summary_out,
            file.path(dir_out_cur_p, "df_summary.csv"),
            quote = F,
            row.names = F)
  
  
  # End of loop on the rows of todo_df
}
