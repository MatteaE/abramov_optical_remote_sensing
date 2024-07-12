# This script mosaics HRS DEMs where the area of interest is
# split across two partially overlapping stereo pairs.
# For Abramov glacier: 2003-08-27, 2005-07-03, 2005-07-29, 2005-09-19.

# Algorithm:
# Crop each half to the Abramov region extent (700000, 731000, 4377000, 4403000)
# Individually co-register to NASADEM (using NK+Tilt; unstable terrain: RGI7.0 + 200 m buffer + NASADEM slopes > 30° + all cells with abs(dh) > 4*NMAD(dh), dh being the pre-coregistration DEMdiff)
# Mosaic the two halves with simple arithmetic average. 
# The results are directly written to folder dem_v4/

library(terra)

# This function takes a DEMdiff SpatRaster and returns
# a SpatVector mask with polygons covering all the
# DEMdiff cells whose absolute value is above 4*NMAD.
func_mad_mask <- function(demdiff) {
  
  mad_cur <- mad(demdiff, na.rm = T) # This is actually the NMAD thanks to default argument "constant = 1.4826"
  dd_mask <- abs(demdiff) > 4*mad_cur
  dd_mask2 <- subst(dd_mask, 0, NA)
  dd_maskv <- as.polygons(dd_mask2)
  return(dd_maskv)
  
}


target_dates <- c("2003-08-27", "2005-07-03", "2005-07-29", "2005-09-19")
dir_data <- "../dem_v3/"
dir_out <- "../dem_v4/"
ext_out <- ext(700000, 731000, 4377000, 4403000)+1000
nasadem_p <- "/PATH/TO/DEM/FOR/COREGISTRATION/"
mask_rgi_slope_p <- "/PATH/TO/MASK/WITH/OUTLINES/AND/STEEP/SLOPES/VECTOR/FILE/"

tmp_d <- "./tmp"
dir.create(tmp_d)
mask_rgi_slope <- vect(mask_rgi_slope_p)


# Prepare NASADEM reference.
rast_blueprint_proc <- extend(crop(rast("/PATH/TO/REFERENCE/GEOTIFF/FOR/GRID/SPECIFICATION/"), ext_out), ext_out)
nasadem_4326 <- rast(nasadem_p)
nasadem_proc <- project(nasadem_4326, rast_blueprint_proc, method = "bilinear")
nasadem_proc_p <- file.path(tmp_d, "nasadem.tif")
writeRaster(nasadem_proc, nasadem_proc_p, overwrite = T)

t_n <- length(target_dates)
for (t_id in 1:t_n) {
  
  cat("Processing", t_id, "/", t_n, "\n")
  
  fns <- paste0(target_dates[t_id], c("a", "b"), "_DEM_30m-adj_coreg.tif")
  
  d1 <- rast(file.path(dir_data, fns[1]))
  d1c <- extend(crop(d1, ext_out), ext_out)
  writeRaster(d1c, file.path(tmp_d, fns[1]), overwrite = T)
  
  d2 <- rast(file.path(dir_data, fns[2]))
  d2c <- extend(crop(d2, ext_out), ext_out)
  writeRaster(d2c, file.path(tmp_d, fns[2]), overwrite = T)
  
  # Prepare unstable mask for first half.
  d1_mask_dh <- func_mad_mask(d1c - nasadem_proc)
  d1_mask <- rbind(d1_mask_dh, mask_rgi_slope)
  writeVector(d1_mask, file.path(tmp_d, "d1_mask.gpkg"))
  
  # Coregister first half.
  cat("Coregistering first half...\n")
  d1_out_p <- file.path(tmp_d, "out1.tif")
  cmd <- "python"
  args <- c("coreg_without_crop.py",
            nasadem_proc_p,
            file.path(tmp_d, fns[1]),
            file.path(tmp_d, "d1_mask.gpkg"),
            d1_out_p)
  system2(cmd, args)
  
  # Prepare unstable mask for second half.
  d2_mask_dh <- func_mad_mask(d2c - nasadem_proc)
  d2_mask <- rbind(d2_mask_dh, mask_rgi_slope)
  writeVector(d2_mask, file.path(tmp_d, "d2_mask.gpkg"))
  
  # Coregister second half.
  cat("Coregistering second half...\n")
  d2_out_p <- file.path(tmp_d, "out2.tif")
  cmd <- "python"
  args <- c("coreg_without_crop.py",
            nasadem_proc_p,
            file.path(tmp_d, fns[2]),
            file.path(tmp_d, "d2_mask.gpkg"),
            d2_out_p)
  system2(cmd, args)
  
  d1cor <- rast(d1_out_p)
  d2cor <- rast(d2_out_p)
  mos <- mosaic(d1cor, d2cor)
  diff_new <- d2cor - d1cor
  writeRaster(diff_new, file.path(tmp_d, paste0("diff_", target_dates[t_id], ".tif")))
  
  writeRaster(mos, file.path(dir_out, paste0(target_dates[t_id], "_DEM_30m-adj_coreg.tif")))
  
  
  file.remove(file.path(tmp_d, "d1_mask.gpkg"),
              file.path(tmp_d, "d2_mask.gpkg"),
              d1_out_p,
              d2_out_p)
  
}
