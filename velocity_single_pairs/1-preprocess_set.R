# This script applies crop and CLAHE (optionally also orientation filter)
# to all IRS-1C/D and SPOT scenes within a directory.
# Filters are performed by MATLAB (code derived from GIV: Van Wyk de Vries and Wickert, 2021).
# called in batch mode on scripts which are generated on the fly.
library(terra)
wd <- "./"
setwd(wd)
source("1a-preprocess_set_matlab_code.R")
source("1b-preprocess_set_matlab_code_only_clahe.R")

apply_matlab_filter <- TRUE # If true, this also applies CLAHE and orientation filter.
matlab_filter_which <- 2 # 1 = orientation + clahe, 2 = only clahe.

ext_sel <- ext(713500, 724500, 4384500, 4394500)

path_input  <- "1-scenes_all"               # Here put the initial .tif orthophotos (inconsistent extents are accepted)
path_output <- "2-scenes_all_cropped_clahe" # Here the program puts the cropped and filtered orthophotos

dir.create(path_output)

path_matlab_bin     <- "PATH/TO/bin/matlab"
path_matlab_script  <- file.path(wd, "matlab_script_cur.m")

lf <- list.files(path_input, pattern = "\\.tif$")
lfn <- length(lf)

for (f_id in 1:lfn) {
  
  f_cur <- lf[f_id]
  
  cat("Processing", f_cur, "--", f_id, "/", lfn, "\n")
  
  p_cur <- file.path(path_input, f_cur)
  
  p_out <- file.path(path_output, f_cur)
  
  # Process only if not already processed
  # (in case we have added files to the selection).
  if (!file.exists(p_out)) {
    
    r_cur <- rast(p_cur)
    
    r_crop <- crop(r_cur, ext_sel)
    
    writeRaster(r_crop, p_out)
    
    if (apply_matlab_filter) {
      
      # Prepare paths for MATLAB script.
      matlab_filter_script_defs <- c(
        paste0("fn='", f_cur, "';"),
        paste0("datafolder='", path_output, "';"),
        paste0("outfolder='", path_output, "';")
      )
      
      # Assemble MATLAB script.
      
      if (matlab_filter_which == 1) {
        matlab_filter_script <- c(matlab_filter_script_defs,
                                  matlab_filter_script_code)
      } else {
        matlab_filter_script <- c(matlab_filter_script_defs,
                                  matlab_filter_script_code_clahe)
      }
      
      # Write MATLAB script.
      writeLines(matlab_filter_script,
                 path_matlab_script)
      
      # Call MATLAB on script, to apply orientation filter and CLAHE.
      system2(path_matlab_bin,
              args = c("-batch",
                       paste0("\"run('", path_matlab_script, "');exit;\"")))
      
      file.remove(path_matlab_script)
      
    }
    
  }
  
}
