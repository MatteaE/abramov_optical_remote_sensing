# Here we load Cosi-CORR results, we compute velocities
# (depending on time separation) and we combine/stack them.

# WARNING: this script must be run with a relatively new version of terra-package,
# because there used to be a bug in the terra::median() function in the CRAN version as of 2023.01.11.


wd <- "./"
setwd(wd)

library(ggplot2)
library(matrixStats)
library(metR)
library(terra)
library(tidyterra)
library(tools)
library(viridis)
source("3b-func_focal_filter_iterate.R")



# Set filtering thresholds.
filter_snr_thresh             <- 0.97 # Minimum SNR to keep a retrieval. From Nanni et al.
filter_vel_mag_thresh         <- 200  # Maximum velocity magnitude [m/yr] to keep a retrieval.
filter_vel_clump_min_size     <- 100    # Remove all velocity retrievals which come in clumps smaller than this size [number of cells].
filter_sd_coeff               <- 2.0    # Coefficient for median+SD filter: remove velocities outside median +/- coeff * SD (or IQR or NMAD)
filter_focal_windowsize       <- 11     # Size of the window for the iterative focal median filter [px]. NOTE it depends on the step of the correlations!
filter_focal_iter_max_n       <- 1      # For the iterative focal filter, limit how many iterations before forced stop. 1 = apply just once.

# Set directory where the displacements are stored.
path_correlations <- file.path(wd, "3-correlations", "scenes_all_clahe_w128_step8", "main_annual")  # Main one, to get the main results.

# Set output directory for the filtered velocity magnitudes.
dir_out <- file.path("4-velocities", "main_annual")      # Main one, to get the main results.
dir.create(dir_out, recursive = TRUE)

# Define stable terrain mask.
# It is a geotiff corresponding to the correlation
# grid, with 0 = unstable, 1 = stable).
path_stable_terrain <- "ref/stable_mask_step8.tif"




#### Load Cosi-CORR results, compute velocities ---------------------------------------------------
# We load grids found in the output directory, based on their time separation.
# This allows us to limit the analysis and aggregation only to some pairs,
# without having to recompute the whole thing.
cat("\n\n==== Loading displacements ====\n")

# List files and load all grids.
lf <- list.files(path_correlations, pattern = "\\.tif$")

# Compute number of characters in the original filenames (except for the file extension).
# filename_nchar_noext <- nchar(format(as.Date("2000/01/01"), format = filename_format_noext))
date_format <- "%Y-%m-%d"

date1_all <- as.Date(substr(lf, 1, 10), format = date_format)
date2_all <- as.Date(substr(regmatches(lf, regexpr("VS_[0-9]{4}-[0-9]{2}-[0-9]{2}", lf)), 4, 13), format = date_format)

lp <- file.path(path_correlations, lf)
displ_all <- rast(lp)

nl <- nlyr(displ_all)
displ_maps_n <- nl / 3

# Split bands into EW, NS and SNR.
displ_ew_rast_all <- displ_all[[seq(1,nl,3)]]
displ_ns_rast_all <- displ_all[[seq(2,nl,3)]]
displ_snr_rast_all <- displ_all[[seq(3,nl,3)]]

# Recompute time separation in days, taking it from the name of each Cosi-CORR output file.
datediff_all <- as.integer(date2_all - date1_all)

# Convert displacement into annual velocity.
vel_ew_rast_all <- displ_ew_rast_all * 365.2425 / datediff_all
vel_ns_rast_all <- displ_ns_rast_all * 365.2425 / datediff_all




#### Simple filtering -----------------------------------------------------------------------------
# We compute filtering masks for all grid layers simultaneously,
# masks have either 1 or NA. Then we multiply layer-wise.

# SNR filter
filter_snr_mask <- subst(displ_snr_rast_all > filter_snr_thresh, 0, NA)
vel_ew_rast_all_f1 <- vel_ew_rast_all * filter_snr_mask
vel_ns_rast_all_f1 <- vel_ns_rast_all * filter_snr_mask


# Remove mean (actually median) EW and NS displacements over stable terrain.
# We do this before magnitude masking, since misregistration
# of 10 m over 10 days already amounts to 365 m/yr bias and needs
# to be preserved (not filtered out) in order to subtract it.
# We use (col)Median(s) and not Mean, since we have not yet filtered for magnitude.
stable_terrain <- rast(path_stable_terrain)
stable_terrain_ids <- which(values(stable_terrain) == 1)

vel_ew_mean_stable <- as.numeric(colMedians(as.matrix(vel_ew_rast_all_f1[stable_terrain_ids]), na.rm = T))
vel_ns_mean_stable <- as.numeric(colMedians(as.matrix(vel_ns_rast_all_f1[stable_terrain_ids]), na.rm = T))

vel_ew_mean_stable[is.na(vel_ew_mean_stable)] <- 0.0
vel_ns_mean_stable[is.na(vel_ns_mean_stable)] <- 0.0

vel_ew_rast_all_f1 <- vel_ew_rast_all_f1 - vel_ew_mean_stable
vel_ns_rast_all_f1 <- vel_ns_rast_all_f1 - vel_ns_mean_stable


# Velocity magnitude filter.
vel_mag_f1 <- sqrt(vel_ew_rast_all_f1^2 + vel_ns_rast_all_f1^2)
filter_vel_mag_mask <- subst(vel_mag_f1 < filter_vel_mag_thresh, 0, NA)
vel_ew_rast_all_f2 <- vel_ew_rast_all_f1 * filter_vel_mag_mask
vel_ns_rast_all_f2 <- vel_ns_rast_all_f1 * filter_vel_mag_mask

vel_ew_rast_all_f2_t <- terra::trim(vel_ew_rast_all_f2)
vel_ns_rast_all_f2_t <- terra::trim(vel_ns_rast_all_f2)



# Filtering of small specks.
# We have to use a for loop because for some reason
# terra::zonal() produces a wrong result on a SpatRaster
# with many layers.
vel_ew_rast_all_f3 <- vel_ew_rast_all_f2_t
vel_ns_rast_all_f3 <- vel_ns_rast_all_f2_t
vel_mag_rast_all_f3 <- sqrt(vel_ew_rast_all_f3*vel_ew_rast_all_f3 + vel_ns_rast_all_f3*vel_ns_rast_all_f3)
for (l_id in 1:length(lf)) {
  cat(l_id, "\n")
  vel_patches <- patches(vel_mag_rast_all_f3[[l_id]],
                         directions = 8,
                         zeroAsNA = FALSE,
                         allowGaps = TRUE)
  vel_zonal <- zonal(setValues(vel_patches, 1), vel_patches, sum, as.raster=TRUE)
  vel_ew_rast_all_f3[[l_id]] <- ifel(vel_zonal < filter_vel_clump_min_size, NA, vel_ew_rast_all_f3[[l_id]])
  vel_ns_rast_all_f3[[l_id]] <- ifel(vel_zonal < filter_vel_clump_min_size, NA, vel_ns_rast_all_f3[[l_id]])
}




# Component-wise, (optionally) iterative focal median filtering of outliers, using focal with median +/- N * SD.
vel_ew_ns_rast_all_f4_l <- func_focal_filter_iterate(vel_ew_rast_all_f3,
                                                     vel_ns_rast_all_f3,
                                                     filter_focal_windowsize,
                                                     filter_sd_coeff,
                                                     filter_focal_iter_max_n)


vel_ew_rast_all_f4 <- vel_ew_ns_rast_all_f4_l[[1]]
vel_ns_rast_all_f4 <- vel_ew_ns_rast_all_f4_l[[2]]

vel_mag_rast_all_f4 <- sqrt(vel_ew_rast_all_f4*vel_ew_rast_all_f4 + vel_ns_rast_all_f4*vel_ns_rast_all_f4)




# Write all velocities. We then extract profiles in the next script.
for (l_id in 1:length(lf)) {
  
  writeRaster(vel_mag_rast_all_f4[[l_id]],      # Main one, to get the main results.
              file.path(dir_out, paste0(format(date1_all[l_id], "%Y-%m-%d_"),
                                        format(date2_all[l_id], "%Y-%m-%d.tif"))),
              overwrite = TRUE)
  
}
