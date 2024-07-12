# This script runs the first processing step of CosiCorr orthorectification for all
# SPOT scenes within the selected target_dir.
# That is, we call tie point detection and selection with several different
# parameter combinations (result is quite dependent on having a good set of GCPs).
# We have two variable parameters: ipfind npoints and ipmatch inlier threshold.
# Then we call CosiCorr RSM generation and refinement, saving the results for each
# combination of parameters.
# Then we will have a second script going through the results and finding the best ones,
# which will be used to actually call the orthorectification.
# When there are both MS and PAN scenes for a certain date, we process both, because
# the PAN scene has higher resolution, but the MS scene (band 1) is sometimes less saturated.
# Then we add a _b1 suffix to the orthorectified MS scene.
# NOTE: this script is hardcoded to work with EPSG:32642 (UTM 42N) - edit as needed.

#### Set paths, modules and environment -----------------------------------------------------------
library(ggplot2)
library(matrixStats)
library(terra)
library(tidyterra)
library(tools)
basepath <- "./"
work_dir   <- normalizePath(file.path(basepath))     # Contains the python orthorectification script and will hold the temporary files during processing
setwd(work_dir)
source("func_process_matches.R")

work_proc_dir <- file.path(work_dir, "proc/")
unlink(work_proc_dir, recursive = TRUE)
dir.create(work_proc_dir)
target_dir <- normalizePath(file.path(basepath, "data/")) # Directory to be processed, it must contain subfolders ms/ and pan/.


env_path <- Sys.getenv("PATH")
Sys.setenv("PATH" = paste0(env_path, ":/PATH/TO/BIN/FOLDER/OF/AMES/STEREO/PIPELINE"))
env_ld_path <- Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv("LD_LIBRARY_PATH" = paste0(env_ld_path, ":/PATH/TO/FOLDER/geoCosiCorr3D/lib"))


python_cmd <- "python3.10"


# Reference data: unstable mask (to filter GCPs), ortho and DEM.
tps_mask_p <- file.path(work_dir, "ref/rgi70.gpkg")
ref_ortho_p <- file.path(work_dir, "ref/ref_ortho.tif")
ref_dem_p <- file.path(work_dir, "ref/ref_dem_egm2008.tif")

ref_ortho_r <- rast(ref_ortho_p)
ref_dem_r <- rast(ref_dem_p)


#### Set parameters for TPs and GCPs --------------------------------------------------------------
ipfind_npoints <- c(1e5, 2.5e5, 5e5, 1e6, 2.5e6, 5e6, 10e6)
ipfind_npoints_dirs <- sapply(ipfind_npoints, function(x) {ifelse(x < 1e6, paste0(x/1000, "k"), paste0(x/1e6, "M"))}) # These are the names of the directories where the feature points of the reference orthos are stored. Sorted by specified number of points.
ipfind_npoints_n <- length(ipfind_npoints)

ipmatch_thresholds <- c(10, 25, 50, 100)
ipmatch_thresholds_n <- length(ipmatch_thresholds)

npoints_prev_max <- 5000 # We use this to limit the increase of the RANSAC inlier threshold on good scenes. If the previous value got us this many valid matches or more, don't further increase the threshold.

ransac_iterations <- 1e6

min_tps  <- 100
max_gcps <- 60


subdirs1 <- normalizePath(list.dirs(target_dir, recursive = FALSE)) # Full path to ms/ and pan/.
#### Processing loop over all scenes --------------------------------------------------------------
for (subdir1_id in 1:length(subdirs1)) {
  
  
  subdir1 <- subdirs1[subdir1_id]
  
  scene_dirs <- list.dirs(subdir1, recursive = FALSE)
  
  
  for (subdir2_id in 1:length(scene_dirs)) {
    
    
    #### . Processing of a single scene -----------------------------------------------------------
    scene_dir <- scene_dirs[subdir2_id]
    
    cat("\n\n\n\n******************************************************************************************\n")
    cat("********** Processing", basename(scene_dir), " **********")
    cat("\n******************************************************************************************\n\n")
    
    date_cur_str <- substr(basename(scene_dir), 22, 31)
    hour_cur_str <- substr(basename(scene_dir), 33, 40)
    
    sensor <- paste0("Spot", substr(basename(scene_dir), 10, 10))
    file.copy(scene_dir, work_proc_dir, recursive = TRUE)
    
    # Path to scene dir as copied into work dir.
    scene_in_work_dir <- file.path(work_proc_dir, basename(scene_dir))
    
    tif_in_work_dir <- file.path(scene_in_work_dir, "SCENE01", "IMAGERY.TIF")
    
    # Name that we will give to the output orthorectified
    # file (without extension): just the date, with an
    # additional _b1 if it comes from an MS scene.
    outname <- paste0(date_cur_str, "_", hour_cur_str)
    
    #### . . If multi-spectral, extract first band ------------------------------------------------
    if (basename(subdir1) == "ms") {
      
      # Extract first band - for all SPOT MS it is the least likely to be saturated.
      # We do not use terra::rast() because it does not like the non-standard
      # georeferencing of the IMAGERY.TIF grid.
      r <- system2("gdal_translate",
                   c("-b", "1",
                     tif_in_work_dir,
                     file.path(scene_in_work_dir, "SCENE01", "IMAGERY_new.TIF")))
      file.rename(file.path(scene_in_work_dir, "SCENE01", "IMAGERY_new.TIF"),
                  tif_in_work_dir)
      
      outname <- paste0(outname, "_b1")
    }
    
    
    #### . . Find GSD based on sensor and mode ----------------------------------------------------
    if (sensor == "Spot5") {
      if (basename(subdir1) == "ms") {gsd <- 10} else {gsd <- 5}
    } else {
      if (basename(subdir1) == "ms") {gsd <- 20} else {gsd <- 10}
    }
    
    
    #### . . Create directories for processing ----------------------------------------------------
    # Main directory with all processing files.
    scene_procdir <- file.path(work_proc_dir, outname)
    dir.create(scene_procdir, recursive = TRUE)
    
    # Here we will put the ipfind feature points of the raw image.
    ipfind_raw_dir <- file.path(scene_procdir, "raw_ipfind")
    unlink(ipfind_raw_dir, recursive = TRUE)
    dir.create(ipfind_raw_dir)
    invisible(sapply(file.path(ipfind_raw_dir, ipfind_npoints_dirs), dir.create))
    
    # Here we will put the ipfind feature points of the cropped ref image.
    ipfind_ref_dir <- file.path(scene_procdir, "ref_ipfind")
    unlink(ipfind_ref_dir, recursive = TRUE)
    dir.create(ipfind_ref_dir)
    invisible(sapply(file.path(ipfind_ref_dir, ipfind_npoints_dirs), dir.create))
    
    # Here we will put the ipmatch results.
    ipmatch_reldir <- file.path(basename(work_proc_dir), outname, "ipmatch/")
    ipmatch_dir <- normalizePath(ipmatch_reldir)
    unlink(ipmatch_dir, recursive = TRUE)
    dir.create(ipmatch_dir)
    
    
    
    #### . . Extract crops of reference ortho and DEM, to speed up processing. --------------------
    # NOTE: the terra::ext() reported by IMAGERY.TIF is wrong due to rotated coordinates,
    #       so we first run rectify() - then the extent is correct.
    # So we use terra::ext(), project the coordinates to UTM, add 5 km of padding on all sides and round to multiple of 10.
    raw_r <- rast(tif_in_work_dir)
    raw_ext <- ext(rectify(raw_r, method = "near"))
    raw_ext_32642 <- project(raw_ext, from = "EPSG:4326", to = "EPSG:32642")
    ext_crop_padding <- 5000
    ref_ext_crop <- ext(c(round(raw_ext_32642[1]/10)*10 - ext_crop_padding,
                          round(raw_ext_32642[2]/10)*10 + ext_crop_padding,
                          round(raw_ext_32642[3]/10)*10 - ext_crop_padding,
                          round(raw_ext_32642[4]/10)*10 + ext_crop_padding))
    ref_ortho_crop_r <- crop(ref_ortho_r, ref_ext_crop)
    if (gsd == 5) {
      ref_ortho_crop_res_r <- disagg(ref_ortho_crop_r, fact = 2, method = "bilinear")
    } else if (gsd == 20) {
      ref_ortho_crop_res_r <- aggregate(ref_ortho_crop_r, fact = 2, fun = "mean")
    } else if (gsd == 10) {
      ref_ortho_crop_res_r <- ref_ortho_crop_r
    } else {
      stop("GSD is ", gsd, " -- something has gone wrong.")
    }
    
    # Rescale the reference to 0-255, with histogram equalization.
    # Equalization helps finding way more keypoints, also appears to help RSM refinement to converge faster)
    ref_ortho_crop_res_stretch_r <- stretch(ref_ortho_crop_res_r, histeq = TRUE) * 255
    
    # Write reference ortho and DEM crops.
    ref_dem_crop_r <- crop(ref_dem_r, ref_ext_crop)
    ref_ortho_final_p <- file.path(scene_procdir, "ref_ortho.tif")
    ref_dem_crop_p <- file.path(scene_procdir, "ref_dem.tif")
    writeRaster(ref_ortho_crop_res_stretch_r, ref_ortho_final_p, overwrite = TRUE)
    writeRaster(ref_dem_crop_r, ref_dem_crop_p, overwrite = TRUE)
    
    
    
    #### . . Iterate over the ipfind npoints ------------------------------------------------------
    # Detect points on the ref and raw images, with the same parameters.
    for (npoints_iter_id in 1:ipfind_npoints_n) {
      
      ipfind_npoints_cur     <- ipfind_npoints[npoints_iter_id]
      ipfind_npoints_str_cur <- ipfind_npoints_dirs[npoints_iter_id]
      
      #### . . . Call ipfind on ref image ---------------------------------------------------------
      ipfind_ref_vwip_dir    <- file.path(ipfind_ref_dir, ipfind_npoints_str_cur)
      cmd <- "ipfind"
      args <- c("--ip-per-image", format(ipfind_npoints_cur, scientific = FALSE),
                "--normalize",
                "--output-folder", ipfind_ref_vwip_dir,
                ref_ortho_final_p)
      system2(cmd, args)
      
      
      #### . . . Call ipfind on raw image ---------------------------------------------------------
      ipfind_raw_vwip_dir    <- file.path(ipfind_raw_dir, ipfind_npoints_str_cur)
      cmd <- "ipfind"
      args <- c("--ip-per-image", format(ipfind_npoints_cur, scientific = FALSE),
                "--normalize",
                "--output-folder", ipfind_raw_vwip_dir,
                tif_in_work_dir)
      system2(cmd, args)
      
      
      ipfind_ref_vwip_p <- file.path(ipfind_ref_vwip_dir, "ref_ortho.vwip")
      ipfind_raw_vwip_p <- file.path(ipfind_raw_vwip_dir, "IMAGERY.vwip")
      
      
      # Now look for point matches with different parameter values.
      # We use match_prev_npoints to avoid raising the threshold too much
      # when a strict threshold already provides plenty of matches
      # (good-quality scenes with plenty of overlap).
      # Else the algorithm to select GCPs has to deal
      # with too many matches (e.g. when computing inter-points
      # distances).
      match_npoints_prev <- 0
      
      #### . . . Iterate over the RANSAC threshold values -----------------------------------------
      for (thresholds_iter_id in 1:ipmatch_thresholds_n) {
        
        if (match_npoints_prev < npoints_prev_max) {
          threshold_cur <- ipmatch_thresholds[thresholds_iter_id]
          
          #### . . . . Call ipmatch ---------------------------------------------------------------
          cmd <- "ipmatch"
          args <- c("-o", paste0(ipmatch_reldir, "match_", ipfind_npoints_str_cur, "_", threshold_cur),
                    "--threads", "8",
                    "--inlier-threshold", threshold_cur,
                    "--ransac-iterations", format(ransac_iterations, scientific = FALSE),
                    ref_ortho_final_p,
                    tif_in_work_dir,
                    ipfind_ref_vwip_p,
                    ipfind_raw_vwip_p)
          system2(cmd, args)
          
          match_file_orig_p <- file.path(ipmatch_dir, paste0("match_", ipfind_npoints_str_cur, "_", threshold_cur, "-", basename(paste0(file_path_sans_ext(ref_ortho_final_p), "__IMAGERY.match"))))
          
          # Did ipmatch produce a match file?
          if (file.exists(match_file_orig_p)) {
            
            #### . . . . Convert match file to text -----------------------------------------------
            match_file_text_p <- file.path(ipmatch_dir, paste0("match_", ipfind_npoints_str_cur, "_", threshold_cur, ".txt"))
            cmd <- python_cmd
            args <- c("proc_match_file.py",
                      match_file_orig_p,
                      match_file_text_p)
            system2(cmd, args, stdout = NULL)
            
            matches_raw <- readLines(match_file_text_p)
            matches_n <- length(matches_raw)
            
            match_npoints_prev <- matches_n # Save this to see whether we can exit the inlier-threshold loop at the next iteration.
            
            # Does the match file have enough points?
            if (matches_n >= min_tps) {
              
              proc_subdir <- file.path(scene_procdir, "proc_versions", paste0(ipfind_npoints_str_cur, "_", threshold_cur))
              dir.create(proc_subdir, recursive = TRUE)
              
              #### . . . . Select GCPs for Cosi-Corr optimization ---------------------------------
              cat("Selecting up to", max_gcps, "points for the Cosi-Corr optimization...\n")
              tps_proc_list <- func_process_matches(ref_ortho_final_p,
                                                    tif_in_work_dir,
                                                    match_file_text_p,
                                                    tps_mask_p,
                                                    gsd,
                                                    max_gcps)
              
              tps_out_p <- file.path(proc_subdir, "tps_proc.txt")
              write.table(tps_proc_list$tps_out,
                          file = tps_out_p,
                          quote = F,
                          col.names = F,
                          row.names = F,
                          sep = "\t")
              cat("Producing GCP plots...\n")
              ggsave(filename = file.path(proc_subdir, "ref.jpg"), plot = tps_proc_list$pl_ref, width = 20, height = 20 * nrow(ref_ortho_crop_res_r) / ncol(ref_ortho_crop_res_r))
              ggsave(filename = file.path(proc_subdir, "raw.jpg"), plot = tps_proc_list$pl_raw, width = 20, height = 20 * nrow(raw_r) / ncol(raw_r))
              
              
              
              
              #### . . . . Run Cosi-Corr optimization ---------------------------------------------
              cat("Running Cosi-CORR RSM refinement...\n")
              cmd <- python_cmd
              args <- c("cosicorr_refinement.py",
                        sensor,
                        gsd,
                        proc_subdir,
                        file.path(scene_in_work_dir, "SCENE01/"),
                        ref_ortho_final_p,
                        ref_dem_crop_p,
                        tps_out_p)
              system2(cmd, args)
              
              # Cleanup: remove the GCP patches.
              # unlink(file.path(proc_subdir, "RSM_refinement", "RSM_gcp_patches"), recursive = TRUE)
              
            } # End do processing if ipmatch gave at least min_tps points.
            
            # After processing, remove temporary match files.
            # file.remove(match_file_orig_p, match_file_text_p)
            
          } # End did ipmatch produce a file?
          
          # End proceed only if previous value of inlier threshold produced not more than npoints_prev_max points (else we are getting way too many points).
        } else {
          break # Else, move on to the next loop over the npoints (restart at lowest value of inlier threshold).
        }
      } # End iteration on RANSAC inlier threshold values in ipmatch.
      
      
    } # End iteration on npoints values in ipfind.
    
    # Cleanup.
    unlink(scene_in_work_dir, recursive = TRUE)
    # Optionally remove the ipfind and ipmatch results.
    # unlink(ipfind_raw_dir, recursive = TRUE)
    # unlink(ipfind_ref_dir, recursive = TRUE)
    # unlink(ipmatch_dir, recursive = TRUE)
    
    
  } # End iteration on the scenes within pan and within ms.
} # End iteration on the subdirs (pan and ms).
