# This script calls MATLAB to apply orientation filter and CLAHE
# to all Sentinel-2 scenes of all sets.
# MATLAB code derived from GIV (Van Wyk de Vries and Wickert, 2021).
# MATLAB script is generated on the fly, loop over all scenes in a set
# is done in MATLAB to avoid restarting MATLAB every time.
# It is much faster.
wd <- "./"
setwd(wd)
source("1a-preprocess_set_matlab_code.R")


# Overview:
# . Loop on orbital path (A and B)
# . . Loop on band (1,2,3,4, corresponding to R, G, B and NIR)
# . . . MATLAB script called on all scenes of the folder of a given band within a given year


# output dir structure:
# 2-filtered
#   path<A,B>
#     <2017-2023>
#       b<1,2,3,4>


targetdir_base <- paste0("/PATH/TO/TREE/OF/SELECTED/SCENES/1-sel")
outputdir_base <- paste0("/PATH/TO/TREE/OF/PREPROCESSED/SCENES/2-filtered")

years_target <- 2017:2023
years_n <- length(years_target)

orb_path_names <- c("pathA", "pathB")

path_matlab_bin     <- "/PATH/TO/bin/matlab"
path_matlab_script  <- file.path(wd, "matlab_script_cur.m") # This is generated on the fly.

for (orb_path_id in 1:2) {
  
  orb_path_cur <- orb_path_names[orb_path_id]
  
  for (year_id in 1:years_n) {
    
    year_cur <- years_target[year_id]
    
    cat("  Year", year_cur, "\n")
    
    path_input_base  <- file.path(targetdir_base, orb_path_cur, year_cur)
    path_output_base <- file.path(outputdir_base, orb_path_cur, year_cur)
    
    # Create directories: by year, with 4 subdirectories for the bands.
    dir.create(path_output_base, recursive = TRUE)
  
    
    for (band_id in 1:4) {
      
      path_input  <- file.path(path_input_base, paste0("b", band_id))
      path_output <- file.path(path_output_base, paste0("b", band_id))
      dir.create(path_output)
      
      # Find files to be processed
      lf <- list.files(path_input, pattern = "\\.tif$")
      
      # Prepare paths for MATLAB script.
      matlab_filter_script_defs <- c(
        paste0("fn_all=[\"", paste(lf, collapse = "\", \""), "\"];"),
        paste0("datafolder='", path_input, "';"), # We used to have path_output here as we extracted the green band from RGB and put it alreadty in path_output. No longer the case.
        paste0("outfolder='", path_output, "';")
      )
      
      
      # Assemble MATLAB script.
      matlab_filter_script <- c(matlab_filter_script_defs,
                                matlab_filter_script_code)
      
      # Write MATLAB script.
      writeLines(matlab_filter_script,
                 path_matlab_script)
      
      # Call MATLAB on script, to apply orientation filter and CLAHE.
      system2(path_matlab_bin,
              args = c("-batch",
                       paste0("\"run('", path_matlab_script, "');exit;\"")))
      
      file.remove(path_matlab_script)
      
    } # End loop on bands
  } # End loop on years
} # End loop on orbital paths
