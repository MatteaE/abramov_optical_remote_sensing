# This script looks at the output from 1-compute_refine_rsm.R
# and looks for the best RSM for each scene.
# This is the one with the lowest mean signed error (as is
# done by normal CosiCorr iterations).

#### Set paths and modules ------------------------------------------------------------------------
library(terra) # For convHull
basepath <- "./"
work_dir   <- normalizePath(file.path(basepath))     # Contains the python orthorectification script and will hold the temporary files during processing
setwd(work_dir)

target_dir <- file.path(work_dir, "proc")

# Algorithm:
# create a summary df for all scenes, then:
# for each scene dir:
# create a df with one row per subdir (each subdir is one parameter combination), a few columns (including avgErr, RMSE, index of best cosicorr iteration, area of convex hull of GCPs)
# for each subdir:
# load the optimization report
# find the number of iterations
# for each iteration, compute avgErr, RMSE and area covered by the convex hull of GCPs
# find lowest meanErr which still respects a minimum area of the convex hull
# write in the per-scene and global df
scene_dirs <- list.files(target_dir, pattern = "[0-9]{4}-[0-9]{2}-[0-9]{2}(_.+)?$")
scene_dirs_n <- length(scene_dirs)

#### Define output data.frame ---------------------------------------------------------------------
df_overall <- data.frame(scene_name           = scene_dirs,
                         ipfind_n_best        = NA_integer_,
                         ipmatch_thresh_best  = NA_real_,
                         cosicorr_iter_best   = NA_integer_,
                         meanerr_best         = NA_real_,
                         rmse_best            = NA_real_,
                         points_hull_area_km2 = NA_real_)

# Here we will store the indices of the scenes to be redone,
# i.e. those with either no valid solutions (i.e. RANSAC always failed) or too high error.
ids_bad <- NULL

#### Iterate on the scenes ------------------------------------------------------------------------
for (scene_id in 1:scene_dirs_n) {
  
  cat("\n\n**** scene_id = ", scene_id, " ****\n")
  
  scene_cur <- scene_dirs[scene_id]
  
  cat("\n\n** scene_cur = ", scene_cur, "\n")
  
  combi_dirs <- list.files(file.path(target_dir, scene_cur, "proc_versions"),
                           pattern = "^[0-9_\\.]{0,3}[kM]{1}_[0-9]{2,}$")
  
  combi_n <- length(combi_dirs)
  
  if (combi_n == 0) {
    message("WARNING: no valid combinations found for ", scene_cur)
    ids_bad <- append(ids_bad, scene_id)
  } else {
    
    df_scene <- data.frame(
      ipfind_n             = rep(NA_integer_, combi_n),
      ipmatch_thresh       = NA_integer_,
      points_hull_area_km2 = NA_real_, # Area (in km2) covered by the convex hull of the GCPs. We use this to exclude solutions with insufficient GCPs coverage (which can have low error but lead to inaccurate orthorectification).
      cosicorr_iter_best   = NA_integer_,
      meanerr              = NA_real_,
      rmse                 = NA_real_)
    
    
    cat("combi_n =", combi_n, "\n")
    
    #### . Iterate on the combinations of ipfind and ipmatch parameters for each scene --------------
    for (combi_id in 1:combi_n) {
      
      cat(". combi_id =", combi_id, "\n")
      
      combi_cur <- combi_dirs[combi_id]
      optim_report_path <- file.path(target_dir, scene_cur, "proc_versions", combi_cur, "RSM_refinement", "tps_proc_GCP_opt.opt_report.csv")
      
      if (!file.exists(optim_report_path)) {
        message("WARNING: file tps_proc_GCP_opt.opt_report.csv does not exist (scene:", scene_cur, "; combi: ", combi_cur, ")")
        next
      }
      
      optim_df <- read.csv(optim_report_path)
      
      iter_nmax <- max(optim_df$nbLoop)
      
      meanerrs   <- rep(NA_real_, iter_nmax+1)
      rmses      <- rep(NA_real_, iter_nmax+1)
      hull_areas <- rep(NA_real_, iter_nmax+1) # We compute the convex hull area of each Cosi-Corr iteration (overkill but most rigorous).
      
      #### . . Iterate on the Cosi-Corr loop iterations for each parameter combination --------------
      for (iter_id in 0:iter_nmax) { # We have to count like Python here! From 0.
        
        p_ids <- which(optim_df$nbLoop == iter_id)
        meanerrs[iter_id + 1]   <- sqrt(mean(optim_df$dxPix[p_ids], na.rm = T)^2 + mean(optim_df$dyPix[p_ids], na.rm = T)^2)
        rmses[iter_id + 1]      <- sqrt(sd(optim_df$dxPix[p_ids], na.rm = T)^2 + sd(optim_df$dyPix[p_ids], na.rm = T)^2)
        hull_areas[iter_id + 1] <- expanse(convHull(vect(cbind(optim_df$Lon[p_ids], optim_df$Lat[p_ids]), crs = "EPSG:4326")), unit = "km") 
        
      }
      
      #### . . Find best Cosi-Corr loop iteration (lowest error) ------------------------------------
      iter_best_id                            <- which.min(meanerrs) # iter_best_id starts from 1. meanerr is still a sqrt(dx^2+dy^2), so only positive, so we can use which.min.
      df_scene$ipfind_n[combi_id]             <- regmatches(combi_cur, regexpr(pattern = "^[0-9_\\.]{0,3}[kM]{1}", combi_cur))
      df_scene$ipmatch_thresh[combi_id]       <- regmatches(combi_cur, regexpr(pattern = "[0-9]{2,}$", combi_cur))
      df_scene$cosicorr_iter_best[combi_id]   <- iter_best_id - 1 # Preserve Python indexing.
      df_scene$meanerr[combi_id]              <- round(meanerrs[iter_best_id], digits = 6)
      df_scene$rmse[combi_id]                 <- round(rmses[iter_best_id], digits = 6)
      df_scene$points_hull_area_km2[combi_id] <- hull_areas[iter_best_id]
      
      sink(file.path(target_dir, scene_cur, "scene_refinement_summary.csv"))
      print(format(df_scene, digits = 6, justify = "right"))
      sink()
    }
    
    
    
    
    #### . Find best parameter combination ----------------------------------------------------------
    # This is the combination with the lowest error,
    # among those whose convex hull of GCPs covers at least X % of the scene.
    # We start with X = 80, then lower it to 70, 60, 50, 40, 30, 20, 10 until a combination is found.
    # All SPOT scenes cover 3600 km2 (60x60).
    # We do this only if we had at least one valid combination.
    if (any(!is.na(df_scene$meanerr))) {
      scene_area <- 3600
      gcp_hull_thresh <- 0.8
      combi_suitable_ids <- which((df_scene$points_hull_area_km2 / scene_area) >= gcp_hull_thresh)
      while (length(combi_suitable_ids) == 0) {
        gcp_hull_thresh <- gcp_hull_thresh - 0.1
        combi_suitable_ids <- which((df_scene$points_hull_area_km2 / scene_area) >= gcp_hull_thresh)
      }
      
      df_scene_suitable_combis <- df_scene[combi_suitable_ids,]
      
      scene_combi_best_id                       <- which.min(df_scene_suitable_combis$meanerr)
      df_overall$scene_name[scene_id]           <- scene_cur
      df_overall$ipfind_n_best[scene_id]        <- df_scene_suitable_combis$ipfind_n[scene_combi_best_id]
      df_overall$ipmatch_thresh_best[scene_id]  <- df_scene_suitable_combis$ipmatch_thresh[scene_combi_best_id]
      df_overall$cosicorr_iter_best[scene_id]   <- df_scene_suitable_combis$cosicorr_iter_best[scene_combi_best_id]
      df_overall$meanerr_best[scene_id]         <- df_scene_suitable_combis$meanerr[scene_combi_best_id]
      df_overall$rmse_best[scene_id]            <- df_scene_suitable_combis$rmse[scene_combi_best_id]
      df_overall$points_hull_area_km2[scene_id] <- df_scene_suitable_combis$points_hull_area_km2[scene_combi_best_id]
    }
  }
  
}


sink(file.path(target_dir, "scenes_all_refinement_summary.csv"))
print(format(df_overall, digits = 6, justify = "right"))
sink()

ids_bad <- unique(append(ids_bad, which(is.na(df_overall$meanerr_best) | (df_overall$meanerr_best > 0.1) | (df_overall$rmse_best > 0.5))))
if (length(ids_bad) > 0) {
  message("found ", length(ids_bad), " scenes with failed RANSAC or high coregistration error:")
  cat(paste0(paste0(df_overall$scene_name[ids_bad], " -- ",
                    df_overall$ipfind_n_best[ids_bad], " -- ",
                    df_overall$ipmatch_thresh_best[ids_bad], " -- ",
                    round(df_overall$meanerr_best[ids_bad], 3), " -- ",
                    round(df_overall$rmse_best[ids_bad], 3)), collapse = "\n"))
}
