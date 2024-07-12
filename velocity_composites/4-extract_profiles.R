# This script extracts longitudinal velocity profiles from the annual composites.

wd <- "./"
setwd(wd)

library(cowplot)
library(ggplot2)
library(ggspatial)
library(terra)
library(tidyterra)
library(reshape2)
library(RStoolbox)
library(sf)
source("4a-func_plot.R")


years <- 2018:2023 # Years that we want to examine. 2017 means 2016-17.

years_n <- length(years)

# Background for profile map.
bg <- rast("/PATH/TO/BACKGROUND/GEOTIFF")

path_annual_mosaics_base <- "PATH/TO/TREE/OF/VELOCITY/MOSAICS"


#### Define profiles for velocity extraction ------------------------------------------------------
points_n <- 1000
# Left branch longitudinal velocity profile: from 3-point spline.
dx_lb <- -200 # To easily move the points around (used only 2 lines below).
dy_lb <- 0
lbranch_pts <- rbind(c(717150+dx_lb, 4388000+dy_lb),
                     c(718500+dx_lb, 4388800+dy_lb),
                     c(719625+dx_lb, 4389500+dy_lb),
                     c(720080+dx_lb, 4390850+dy_lb),
                     c(720130+dx_lb, 4391720+dy_lb),
                     c(720303+dx_lb, 4392154+dy_lb))
lbranch_line_coords <- cbind(spline(x = lbranch_pts[,1], n = points_n, method = "natural")$y,
                             spline(x = lbranch_pts[,2], n = points_n, method = "natural")$y)


# Right branch (roughly approximated) flowline.
dx_rb <- 0
dy_rb <- 0
rbranch_pts <- rbind(c(720937+dx_rb, 4386925+dy_rb),
                     c(720000+dx_rb, 4388500+dy_rb),
                     c(719843+dx_rb, 4389493+dy_rb),
                     c(720043+dx_rb, 4390493+dy_rb))
rbranch_line_coords <- cbind(spline(x = rbranch_pts[,1], n = points_n)$y,
                             spline(x = rbranch_pts[,2], n = points_n)$y)


# Prepare profiles for maps.
prof_line_lbranch <- st_sfc(st_linestring(lbranch_line_coords, dim = "XY"), crs = "EPSG:32642")
prof_line_rbranch <- st_sfc(st_linestring(rbranch_line_coords, dim = "XY"), crs = "EPSG:32642")

# Compute along-profile distances.
lbranch_st_points <- st_cast(prof_line_lbranch, "POINT")
lbranch_dists <- c(0.0,cumsum(st_distance(lbranch_st_points[-1],lbranch_st_points[-length(lbranch_st_points)], by_element=TRUE)))
rbranch_st_points <- st_cast(prof_line_rbranch, "POINT")
rbranch_dists <- c(0.0,cumsum(st_distance(rbranch_st_points[-1],rbranch_st_points[-length(rbranch_st_points)], by_element=TRUE)))

# Write profiles as GPKG.
prof_out <- rbind(st_as_sf(prof_line_lbranch),
                  st_as_sf(prof_line_rbranch))
st_write(prof_out, "prof_lines.gpkg", append = FALSE, delete_layer = TRUE)

vel_lbranch_t1 <- data.frame(p_id = 1:points_n, track = 1, dist_along = lbranch_dists)
for (year_id in 1:years_n) {
  vel_lbranch_t1[paste0("year", year_id)] <- NA_real_
}
vel_rbranch_t1 <- vel_lbranch_t1
vel_rbranch_t1$dist_along <- rbranch_dists
vel_lbranch_t2 <- vel_lbranch_t1
vel_lbranch_t2$track <- 2
vel_rbranch_t2 <- vel_lbranch_t2
vel_rbranch_t2$dist_along <- rbranch_dists


# Load velocity grids.
grids_n <- 2 * years_n
grids <- list()
for (year_id in 1:years_n) {
  year_cur <- years[year_id]
  grids[[(year_id-1)*2 + 1]] <- rast(file.path(path_annual_mosaics_base,
                                               "pathA",
                                               year_cur,
                                               "composite_mag.tif"))
  grids[[year_id*2]] <- rast(file.path(path_annual_mosaics_base,
                                       "pathB",
                                       year_cur,
                                       "composite_mag.tif"))
}



#### . Compute stable-terrain NMAD for uncertainties ----------------------------------------------
stable_ras <- crop(rast("stable_mask_40m.tif"), grids[[1]])
vel_uncertainty <- matrix(NA_real_, ncol = 2, nrow = years_n)

# Extract NMAD for each year.
# First column track 1 or A, second column track 2 or B.
for (year_id in 1:years_n) {
  vel_uncertainty[year_id,1] <- mad(values(mask(grids[[((year_id-1)*2)+1]], stable_ras, maskvalues = 0), mat = FALSE, na.rm = T))
  vel_uncertainty[year_id,2] <- mad(values(mask(grids[[(year_id)*2]], stable_ras, maskvalues = 0), mat = FALSE, na.rm = T))
}


#### . Extract velocities and put into data frames ------------------------------------------------
for (year_id in 1:years_n) {
  vel_lbranch_t1[,year_id+3] <- extract(grids[[((year_id-1)*2)+1]], lbranch_line_coords, method = "bilinear")
  vel_rbranch_t1[,year_id+3] <- extract(grids[[((year_id-1)*2)+1]], rbranch_line_coords, method = "bilinear")
  vel_lbranch_t2[,year_id+3] <- extract(grids[[(year_id)*2]], lbranch_line_coords, method = "bilinear")
  vel_rbranch_t2[,year_id+3] <- extract(grids[[(year_id)*2]], rbranch_line_coords, method = "bilinear")
}


# For plots with independent orbital tracks.
vel_lbranch_t1_melt <- melt(vel_lbranch_t1, id.vars = c("p_id", "track", "dist_along"))
vel_lbranch_t2_melt <- melt(vel_lbranch_t2, id.vars = c("p_id", "track", "dist_along"))

vel_rbranch_t1_melt <- melt(vel_rbranch_t1, id.vars = c("p_id", "track", "dist_along"))
vel_rbranch_t2_melt <- melt(vel_rbranch_t2, id.vars = c("p_id", "track", "dist_along"))



#### Now produce the plots ------------------------------------------------------------------------
# Four plots: we separate by branch and by orbital track.
ext_left <- c(715750, 720750, 4387250, 4392250)
ext_right <- c(716750, 722250, 4386250, 4391750)
map_lbranch <- func_plot_profile_map(prof_line_lbranch,
                                     lbranch_pts,
                                     extent = ext_left,
                                     labels = c("1", "3"))
map_rbranch <- func_plot_profile_map(prof_line_rbranch,
                                     rbranch_pts,
                                     extent = ext_right,
                                     labels = c("4", "5"))


# Add uncertainty to the melted df (example on orbital track 1 only).
vel_lbranch_t1_melt$vel_uncertainty <- NA_real_
vel_rbranch_t1_melt$vel_uncertainty <- NA_real_
for (year_id in 1:years_n) {
  ids_year <- which(vel_lbranch_t1_melt$variable == paste0("year", year_id))
  vel_lbranch_t1_melt$vel_uncertainty[ids_year] <- vel_uncertainty[year_id,1]
  vel_rbranch_t1_melt$vel_uncertainty[ids_year] <- vel_uncertainty[year_id,1]
}




p_lines_s2_lbranch <- func_plot_annual_profiles(vel_lbranch_t1_melt, years, labels = c("1", "3"))
p_lines_s2_rbranch <- func_plot_annual_profiles(vel_rbranch_t1_melt, years, labels = c("4", "5"))

save(p_lines_s2_lbranch, p_lines_s2_rbranch,
     file = "sentinel2_speedup_plots.RData")
