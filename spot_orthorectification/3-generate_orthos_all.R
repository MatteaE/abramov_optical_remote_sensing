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

work_proc_dir   <- file.path(work_dir, "proc")
target_dir      <- normalizePath(file.path(basepath, "data/"))              # Directory to be processed, it must contain subfolders ms/ and pan/.
output_all_dir  <- file.path(work_dir, "proc")                              # We store all versions of each orthophoto in the same proc dir as the rest of processing.
output_best_dir <- normalizePath(file.path(basepath, "data/ortho_best/"))   # Directory to store the best version of all orthophotos.
unlink(output_best_dir, recursive = TRUE)
dir.create(output_best_dir, recursive = TRUE)


env_path <- Sys.getenv("PATH")
Sys.setenv("PATH" = paste0(env_path, ":/PATH/TO/BIN/FOLDER/OF/AMES/STEREO/PIPELINE"))
env_ld_path <- Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv("LD_LIBRARY_PATH" = paste0(env_ld_path, ":/PATH/TO/FOLDER/geoCosiCorr3D/lib"))


python_cmd <- "python3.10"



#### Load table with info on which is the expected best parameter combination for each scene ------
df_overall <- read.table(file.path(work_proc_dir, "scenes_all_refinement_summary.csv"), header = T)


subdirs1 <- normalizePath(list.dirs(target_dir, recursive = FALSE)) # Full path to ms/ and pan/.
#### Processing loop over pan/ms ------------------------------------------------------------------
for (subdir1_id in 1:length(subdirs1)) {
  
  
  subdir1 <- subdirs1[subdir1_id]
  
  scene_dirs <- list.dirs(subdir1, recursive = FALSE)
  
  #### . Processing loop over all scenes (either pan or ms) ---------------------------------------
  for (subdir2_id in 1:length(scene_dirs)) {
    
    
    #### . . Processing of a single scene ---------------------------------------------------------
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
    
    
    # This is the directory where 1-compute_refine_rsm.R did
    # its processing (ipfind, ipmatch, RSM refinement).
    scene_procdir <- file.path(work_proc_dir, outname)
    scene_out_all_dir <- file.path(output_all_dir, outname, "ortho/")
    dir.create(scene_out_all_dir, recursive = TRUE)
    
    scene_combi_df <- read.table(file.path(scene_procdir, "scene_refinement_summary.csv"), header = T)
    
    combi_dirs <- list.files(file.path(scene_procdir, "proc_versions"),
                             pattern = "^[0-9_\\.]{0,3}[kM]{1}_[0-9]{2,}$")
    combi_n <- length(combi_dirs)
    
    #### . . . Processing loop over the parameter combinations for a single scene -----------------
    for (combi_id in 1:combi_n) {
      cat("Parameter combination", combi_id, "out of", combi_n, "--")
      
      combi_cur <- combi_dirs[combi_id]
      
      combi_cur_npoints_str <- regmatches(combi_cur, regexpr("^[0-9_\\.]{0,3}[kM]{1}", combi_cur))
      combi_cur_thresh <- as.numeric(regmatches(combi_cur, regexpr(pattern = "[0-9]{2,}$", combi_cur)))
      
      combi_cur_in_df_id <- which((scene_combi_df$ipfind_n == combi_cur_npoints_str) & (scene_combi_df$ipmatch_thresh == combi_cur_thresh))
      if (length(combi_cur_in_df_id) == 1) {
        
        # Find the CosiCorr iteration that produced the lowest error
        # for the current parameter combination. Counted as in Python (from 0).
        # We use it to run orthorectification with its parameters.
        combi_cur_iter_best_id <- scene_combi_df$cosicorr_iter_best[combi_cur_in_df_id] 
        
        output_path <- file.path(scene_out_all_dir, paste0(outname, "_", combi_cur, ".tif"))
        
        cat(" Running orthorectification...\n")
        cmd <- python_cmd
        args <- c("./orthorectify_spot_single.py",
                  "-o", output_path,
                  sensor,
                  gsd,
                  scene_in_work_dir,
                  file.path(scene_procdir, "proc_versions"),
                  file.path(scene_procdir, "ref_dem.tif"),
                  combi_cur,
                  combi_cur_iter_best_id)
        #### . . . . Call single orthorectification -----------------------------------------------
        system2(cmd,args)
        
      } else {
        stop(paste0("could not find combi_cur within scene_combi_df. combi_cur = ", combi_cur, " -- outname = ", outname))
      }
      
    } # End iterate on the parameter combinations for one scene
    
    
    #### . . . Copy best orthorectification result to output directory ----------------------------
    df_overall_scene_id <- which(df_overall$scene_name == outname)
    scene_best_combi_str <- paste0(df_overall$ipfind_n_best[df_overall_scene_id], "_", df_overall$ipmatch_thresh_best[df_overall_scene_id])
    scene_best_ortho_p <- file.path(scene_out_all_dir, paste0(outname, "_", scene_best_combi_str, ".tif"))
    file.copy(scene_best_ortho_p, output_best_dir)
    
    
    #### . . . Remove raw scene after all its orthorectifications are done ------------------------
    unlink(scene_in_work_dir, recursive = TRUE)
    
    
  } # End iterate on the scene dirs
} # End iterate on the subdir (ms / pan)
