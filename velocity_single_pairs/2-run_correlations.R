# This script takes a directory full of pre-processed orthoimages
# called %Y%m%d.tif and runs Cosi-CORR correlation on the specified pairs.
# We compare each image only with later images.
library(terra)
library(tools)

# Set LD_LIBRARY_PATH
env_ld_path <- Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv("LD_LIBRARY_PATH" = paste0(env_ld_path, ":/PATH/TO/geoCosiCorr3D/lib/"))
wd <- "./"
setwd(wd)


#### Set processing options -----------------------------------------------------------------------

# Set paths.
path_input <- "./2-scenes_all_cropped_clahe"

# Main annual results.
path_correlations <- "./3-correlations/scenes_all_clahe_w128_step8/main_annual"


# Select scene pairs to be processed. First column = first scene of the correlation, second column = second scene of the correlation.
corr_sel <- rbind(
  
  # # Annual
  # c("1996-09-15","1997-08-30"),
  # c("1997-08-30","1998-07-28"),
  # c("1998-07-28","1999-07-16"),
  # c("2000-08-14","2001-08-04"),
  
)


#### Iterate on the scenes, compute displacements with Cosi-CORR ----------------------------------
dir.create(path_correlations)

lf <- list.files(path_input, pattern = "^[0-9]{4}-[0-9]{2}-[0-9]{2}(.)*\\.tif$")

corr_n <- nrow(corr_sel)

corr_lf <- list.files(path_correlations, pattern = "\\.tif$")
for (corr_id in 1:corr_n) {
  
  corr_cur <- corr_sel[corr_id,]
  
  cat("==== Processing", paste0(corr_cur, collapse = " vs "), "====\n")
  f1 <- lf[grep(corr_cur[1], lf)]
  f2 <- lf[grep(corr_cur[2], lf)]  
  
  p1 <- file.path(path_input, f1)
  p2 <- file.path(path_input, f2)
  
  out_exists <- intersect(grep(corr_cur[1], corr_lf), grep(corr_cur[2], corr_lf))
  if (length(out_exists) > 0) {
    cat("Correlation file exists, skipping...\n")
  } else {
    
    
    # Call Cosi-CORR.
    system2("python3.10",
            args = c("cosicorr_displacement.py",
                     p1,
                     p2,
                     path_correlations))
  }
  
} # End iteration on all scenes with Cosi-CORR correlation.
