# This script processes the SPOT HRS DEMs to remove
# areas where the mapprojected orthoimage is saturated.
# These correspond to places where the SPOT scene
# had no data, and all data can only come from the seed DEM.
# NOTE: we always mapproject on the seed Copernicus GLO-90 DEM,
# since that's what we used to compute the stereo DEM in the first place.
library(terra)

Sys.setenv("PATH" = paste0(Sys.getenv("PATH"), "/PATH/TO/STEREO/PIPELINE/bin"))

dem_v1_dirpath <- "../dem_v1/"
dem_v2_dirpath <- "../dem_v2/"
repo_dirpath <- "/PATH/TO/TREE/OF/STEREO/PAIRS/"

dem_seed_path <- "/PATH/TO/SEED/DEM/GEOTIFF"
proj4 <- "'+proj=utm +zone=42 +datum=WGS84 +units=m +no_defs +type=crs'"

dems_lf <- list.files(dem_v1_dirpath, pattern = "\\.tif$")
dems_n <- length(dems_lf)

repo_ld <- list.dirs(repo_dirpath, recursive = F, full.names = F)

for (dem_id in 1:dems_n) {
# dem_id <- 1

  demv1_cur_path <- file.path(dem_v1_dirpath, dems_lf[dem_id])
  demv1_cur <- rast(demv1_cur_path)
  
  cat("==========\n Working on", basename(demv1_cur_path), "...\n==========\n")
  
  date_cur_str <- substr(dems_lf[dem_id], 1, 10)
  date_cur <- as.Date(date_cur_str, format = "%Y-%m-%d")
  
  
  # Now mapproject the front and back images.
  pair_cur_path <- file.path(repo_dirpath, repo_ld[dem_id])
  
  pair_cur_ld <- list.dirs(pair_cur_path, recursive = F, full.names = F)
  front_dirpath <- pair_cur_ld[grep("HRS-1", pair_cur_ld, fixed = TRUE)]
  back_dirpath <- pair_cur_ld[grep("HRS-2", pair_cur_ld, fixed = TRUE)]
  
  
  front_dimpath <- file.path(pair_cur_path, front_dirpath, "SCENE01", "METADATA.DIM")
  back_dimpath <- file.path(pair_cur_path, back_dirpath, "SCENE01", "METADATA.DIM")
  
  front_tifpath <- file.path(pair_cur_path, front_dirpath, "SCENE01", "IMAGERY.TIF")
  back_tifpath <- file.path(pair_cur_path, back_dirpath, "SCENE01", "IMAGERY.TIF")
  
  cmd <- "add_spot_rpc"
  args <- c(front_dimpath, "-o", front_dimpath)
  system2(cmd, args)
  
  args <- c(back_dimpath, "-o", back_dimpath)
  system2(cmd, args)
  
  front_out_path <- file.path(pair_cur_path, "front_mapproj.tif")
  cmd <- "mapproject"
  args <- c("-t", "rpc",
            "--t_srs", proj4,
            "--mpp", 10,
            dem_seed_path,
            front_tifpath,
            front_dimpath,
            front_out_path)
  system2(cmd, args)
  
  back_out_path <- file.path(pair_cur_path, "back_mapproj.tif")
  cmd <- "mapproject"
  args <- c("-t", "rpc",
            "--t_srs", proj4,
            "--mpp", 10,
            dem_seed_path,
            back_tifpath,
            back_dimpath,
            back_out_path)
  system2(cmd, args)
  
  # Read in the mapproj image as rast(), vectorize the areas with saturated pixels, use them to mask the DEM, write the updated DEM to v2.
  front_rast <- rast(front_out_path)
  back_rast <- rast(back_out_path)
  
  cat("\n\nComputing saturated masks...\n")
  front_sat_mask <- front_rast >= 255
  back_sat_mask <- back_rast >= 255
  
  bsm_align <- resample(back_sat_mask, front_sat_mask, method = "near")
  
  cat("Computing final saturated mask...\n")
  saturated_mask <- front_sat_mask | bsm_align
  saturated_mask_align <- resample(saturated_mask, demv1_cur, method = "near")
  
  cat("Masking DEM...\n")
  demv2_cur <- mask(demv1_cur, saturated_mask_align, maskvalues = 1)
  
  cat("Writing masked DEM...\n")
  writeRaster(demv2_cur, file.path(dem_v2_dirpath, dems_lf[dem_id]))
  
  # unlink(pair_cur_path, recursive = TRUE)
  
}

