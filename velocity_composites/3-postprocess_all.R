# Here we load Cosi-CORR results, we compute velocities
# (depending on time separation) and we combine/stack them.

wd <- "./"
setwd(wd)

library(ggplot2)
library(matrixStats)
library(metR)
library(terra)
library(tidyterra)
library(tools)
library(viridis)
source("3a-func_postprocess_set.R")
source("3b-func_plot_velocity.R")

# Overview:
# . Loop on orbital path (A and B)
# . . Loop on year (1 to 7, i.e. 2016/17 to 2022/23)
# . . . Loop on band (1,2,3,4, corresponding to R, G, B and NIR)

# Output:
# velocity_mosaics/
# . path<A,B>/
# . . <2017-2023>/
# . . . b<1-4>_<ns,ew,mag,count,frac>.tif
# . . . b<1-4>_plot.pdf
# . . . composite_<ns,ew,mag>.tif


# Set path to a file with rejected dates.
# The dates refer to original Sentinel-2 scenes:
# it is a simple way of filtering out scenes without
# having to redo correlations.
# If a blacklisted date appears in a correlation, the
# correlation is discarded when we create the mosaics.
path_scenes_blacklist <- "scenes_blacklist.txt"      # Blacklist with the dates of all scenes with 10 % or more clouds on Abramov.

# Set filtering thresholds.
filter_snr_thresh     <- 0.97 # Minimum SNR to keep a retrieval. From Nanni et al. preprint.
filter_vel_mag_thresh <- 100  # Maximum velocity magnitude [m/yr] to keep a retrieval.

# Set time specifications to select only some scene pairs.
# These are all the allowed time separations in days.
dt_sel <- c(1:100, 300:430)

# Set directory where the displacements are stored.
path_correlations_base <- "/PATH/TO/TREE/OF/CORRELATIONS"

# Set output directory for velocity mosaics.
path_output_base       <- "PATH/TO/TREE/OF/VELOCITY/MOSAICS"

# Set filename format of the original scenes.
filename_format_noext <- "%Y%m%d"

# Set suffix of Cosi-CORR files.
displ_suffix <- "_frequency_wz_64_step_4.tif"

# Define stable terrain mask.
# It is a geotiff corresponding to the correlation
# grid, with 0 = unstable, 1 = stable).
path_stable_terrain <- "stable_mask_40m.tif"

path_outlines <- "/PATH/TO/RGI/VECTOR/FILE"
outl <- vect(path_outlines)


years_target <- 2017:2023
years_n <- length(years_target)

orb_path_names <- c("pathA", "pathB")

blacklist_dates_str <- readLines(path_scenes_blacklist)


#### Loops start here -----------------------------------------------------------------------------
for (orb_path_id in 1:2) {
  # for (orb_path_id in 1:1) {
  
  orb_path_cur <- orb_path_names[orb_path_id]
  cat("\n==== Orbital path", orb_path_cur, "====\n")
  
  
  for (year_id in 1:years_n) {
    
    year_cur <- years_target[year_id]
    cat("\n . Year", year_cur, "\n")
    
    bands_output_all <- list(ns    = list(),
                             ew    = list(),
                             count = list())
    
    for (band_id in 1:4) {
      
      cat(". . Band", band_id, "\n")
      
      
      path_correlations <- file.path(path_correlations_base, orb_path_cur, year_cur, paste0("b", band_id))
      path_dir_output   <- file.path(path_output_base, orb_path_cur, year_cur)
      
      dir.create(path_dir_output, recursive = TRUE, showWarnings = FALSE)
      
      
      paths_out <-  file.path(path_dir_output, paste0("b", band_id, c("_mag.tif",
                                                                      "_ns.tif",
                                                                      "_ew.tif",
                                                                      "_count.tif",
                                                                      "_frac.tif",
                                                                      "_plot.pdf")))
      
      # If all output files already exist, load them instead of recomputing and writing them.
      if (any(!file.exists(paths_out))) {
        
        band_output_l <- func_postprocess_set(path_correlations,
                                              filename_format_noext,
                                              dt_sel,
                                              blacklist_dates_str)
        
        
        writeRaster(band_output_l$vel_mag,   paths_out[1], overwrite = TRUE)
        writeRaster(band_output_l$vel_ns,    paths_out[2], overwrite = TRUE)
        writeRaster(band_output_l$vel_ew,    paths_out[3], overwrite = TRUE)
        writeRaster(band_output_l$vel_count, paths_out[4], overwrite = TRUE)
        writeRaster(band_output_l$vel_frac,  paths_out[5], overwrite = TRUE)
        ggsave(paths_out[6],
               plot = band_output_l$pl_vel, width = 12, height = 8)
        
      } else {
        
        band_output_l <- list(vel_mag   = rast(paths_out[1]),
                              vel_ns    = rast(paths_out[2]),
                              vel_ew    = rast(paths_out[3]),
                              vel_count = rast(paths_out[4]),
                              vel_frac  = rast(paths_out[5]))
      }
      
      # Save band output for later composite production.
      bands_output_all[["ns"]][[band_id]]    <- band_output_l$vel_ns
      bands_output_all[["ew"]][[band_id]]    <- band_output_l$vel_ew
      bands_output_all[["count"]][[band_id]] <- band_output_l$vel_count
      
    }

    
    # Now compute composite NS and EW of the four bands.
    # Arithmetic average of the four results, weighted
    # (for each pixel) by the scene count of each band.
    # NOTE: we have also tested median, virtually no change vs weighted average.
    cat("\nComputing average of the 4 bands...\n")
    composite_ns  <-  sum(rast(bands_output_all[["ns"]]) * rast(bands_output_all[["count"]])) / sum(rast(bands_output_all[["count"]]))
    composite_ew  <-  sum(rast(bands_output_all[["ew"]]) * rast(bands_output_all[["count"]])) / sum(rast(bands_output_all[["count"]]))
    composite_mag <-  sqrt(composite_ew^2 + composite_ns^2)
    pl_composite <- func_plot_velocity(composite_ew,
                                       composite_ns,
                                       composite_mag)
    
    # Plot and save composite grids.
    writeRaster(composite_mag, file.path(path_dir_output, "composite_mag.tif"), overwrite = TRUE)
    writeRaster(composite_ns,  file.path(path_dir_output, "composite_ns.tif"), overwrite = TRUE)
    writeRaster(composite_ew,  file.path(path_dir_output, "composite_ew.tif"), overwrite = TRUE)
    ggsave(file.path(path_dir_output, "composite_plot.pdf"),
           plot = pl_composite, width = 12, height = 8)
    
  }
}
