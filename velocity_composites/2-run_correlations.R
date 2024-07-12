# This script takes a directory full of pre-processed orthoimages
# called %Y%m%d.tif and runs Cosi-CORR correlation on the specified pairs.
# We compare each image only with later images.
library(terra)
library(tools)

# Overview:
# . Loop on orbital path (A and B)
# . . Loop on year (1 to 7, i.e. 2016/17 to 2022/23)
# . . . Loop on band (1,2,3,4, corresponding to R, G, B and NIR)


env_ld_path <- Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv("LD_LIBRARY_PATH" = paste0(env_ld_path, ":/PATH/TO/geoCosiCorr3D/lib/"))
wd <- "./"
setwd(wd)


#### Set processing options -----------------------------------------------------------------------

inputdir_base        <- paste0("/PATH/TO/TREE/OF/FILTERED/SCENES/2-filtered")
correlationsdir_base <- paste0("/PATH/TO/TREE/OF/CORRELATIONS")

years_target <- 2017:2023
years_n <- length(years_target)

orb_path_names <- c("pathA", "pathB")

# Set filename format.
filename_format_noext <- "%Y%m%d"
filename_ext <- ".tif"

# Set time specifications to select scene pairs.
# These are all the allowed time separations in days.
dt_sel <- c(1:100, 300:430)


for (orb_path_id in 1:2) {
  
  orb_path_cur <- orb_path_names[orb_path_id]
  
  for (year_id in 1:years_n) {
    
    year_cur <- years_target[year_id]
    
    
    for (band_id in 1:4) {
      
      path_input        <- file.path(inputdir_base, orb_path_cur, year_cur, paste0("b", band_id))
      path_correlations <- file.path(correlationsdir_base, orb_path_cur, year_cur, paste0("b", band_id))
      dir.create(path_correlations, recursive = TRUE, showWarnings = FALSE)
      
      lf <- list.files(path_input, pattern = "\\.tif$")
      filename_format <- paste0(filename_format_noext, "_b", band_id, filename_ext)
      f_dates <- as.Date(lf, filename_format)
      lfn <- length(lf)
      
      # Iterate on all scenes, for each find which correlations we want to do and do them.
      for (f_id in 1:lfn) {
        
        f_cur <- lf[f_id]
        
        cat("==== Processing", f_cur, "====\n")
        
        f_date_cur <- f_dates[f_id]
        f_corr_ids <- which((f_dates - f_date_cur) %in% dt_sel)
        
        # How many correlations do we have to do for the current scene?
        corr_curscene_n <- length(f_corr_ids)
        corr_curscene_todo_n <- corr_curscene_n # This we will update by removing already existing correlations.
        cat("   ", corr_curscene_n, "correlations with the current scene.")
        
        # Do any correlations exist already? Skip them.
        p_cur <- file.path(path_input, f_cur)
        p_corr_cur_all  <- normalizePath(file.path(path_input, lf[f_corr_ids]))
        p_corr_cur_todo <- p_corr_cur_all
        
        paths_output_all <- file.path(path_correlations,
                                      paste0(file_path_sans_ext(basename(p_cur)),
                                             "_VS_",
                                             file_path_sans_ext(basename(p_corr_cur_all)),
                                             "_frequency_wz_64_step_4.tif"))
        
        corr_curscene_existing_ids <- which(file.exists(paths_output_all))
        corr_curscene_existing_n   <- length(corr_curscene_existing_ids)
        
        if (corr_curscene_existing_n > 0) {
          p_corr_cur_todo <- p_corr_cur_all[-corr_curscene_existing_ids]
          corr_curscene_todo_n <- corr_curscene_n - corr_curscene_existing_n
          cat("", corr_curscene_existing_n, "exist already and will be skipped.")
        }
        cat("\n")
        
        
        # Proceed only if there are correlations
        # to do for the current scene,
        # else move on to the next.
        # We will call Python just once for each "first" scene
        # of the correlation.
        if (corr_curscene_todo_n > 0) {
          
          p_corr_listfile_p <- "list_to_correlate.txt"
          
          writeLines(p_corr_cur_todo, p_corr_listfile_p)
          
          # Run all correlations.
          # Just one call to Python, faster than iterating on each correlation.
          cmd <- "python"
          args <- c("cosicorr_displacement.py",
                    p_cur,
                    p_corr_listfile_p,
                    path_correlations)
          system2(cmd, args)
          
          # Cleanup.
          file.remove(p_corr_listfile_p)
          
        
        } # End if (corr_curscene_todo_n  > 0).
      } # End iteration on all scenes with Cosi-CORR correlation.
      
    } # End iteration on bands.
    
  } # End iteration on years.
  
} # End iteration on orbital paths.
