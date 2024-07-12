# This script extracts longitudinal profiles of glacier surface velocity,
# filters them and plots them.
library(cowplot)
library(ggplot2)
library(ggspatial)
library(png)
library(terra)
library(tidyterra)
library(reshape2)
library(RStoolbox)
library(sf)
source("4a-func_process_profile.R")

# To select between different versions.
version_suffix <- ""

# Directory with the velocity maps as processed by 3-postprocess_set.R.
vel_maps_dir <- file.path("4-velocities", "main_annual")                         # Main one, to get the main results.


# Set this to NA to select all velocity maps for plot.
# Else it selects velocity maps from the list.files() of the vel_map_dir.
vel_maps_sel_ids <- NA



# Background for profile map.
bg <- rast("/PATH/TO/BACKGROUND/OF/MAP.tif")

#### Define profiles for velocity extraction ------------------------------------------------------
points_n <- 1000
# Left branch longitudinal velocity profile: from spline.
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


#### Processing starts here -----------------------------------------------------------------------

# Prepare profiles for maps.
prof_line_lbranch <- st_sfc(st_linestring(lbranch_line_coords, dim = "XY"), crs = "EPSG:32642")

prof_lines <- st_as_sf(prof_line_lbranch)
st_write(prof_lines, "prof_lines.gpkg", append = FALSE)

lbranch_points <- st_cast(prof_line_lbranch, "POINT")
lbranch_dists <- c(0.0,cumsum(st_distance(lbranch_points[-1],lbranch_points[-length(lbranch_points)], by_element=TRUE)))


vel_maps <- list.files(vel_maps_dir, pattern = "\\.tif$")

if (!is.na(vel_maps_sel_ids[1])) {
  vel_maps <- vel_maps[vel_maps_sel_ids]
}

vel_maps_p <- file.path(vel_maps_dir, vel_maps)
vel_maps_df <- data.frame(filename = vel_maps,
                          filepath = vel_maps_p,
                          name     = tools::file_path_sans_ext(vel_maps),
                          t1       = as.Date(substr(vel_maps, 1, 10)),
                          t2       = as.Date(substr(vel_maps, 12, 21)))
vel_maps_df$dt <- vel_maps_df$t2 - vel_maps_df$t1


# Load all velocity maps.
grids_n <- length(vel_maps)
grids <- list()
for (grid_id in 1:grids_n) {
  grids[[grid_id]] <- rast(vel_maps_p[grid_id])
}

# Create, initialize and fill data frames with the extracted velocity profiles.
vel_lbranch <- data.frame(p_id = 1:points_n,
                          dist_along = lbranch_dists)
for (grid_id in 1:grids_n) {
  vel_lbranch[vel_maps_df$name[grid_id]] <- NA_real_
}
for (grid_id in 1:grids_n) {
  cat("Extracting from grid", grid_id, "/", grids_n, "\n")
  vel_lbranch[,grid_id+2] <- extract(grids[[grid_id]], lbranch_line_coords, method = "simple")
}


# Select same start for all profiles (8-bit sensors do
# not work in the accumulation area due to saturation).
profile_crop_upper <- 3250
ids_upper <- which(vel_lbranch$dist_along < profile_crop_upper)
vel_lbranch[ids_upper,3:ncol(vel_lbranch)] <- NA_real_



# Process and filter velocity profiles
cat("Pruning and filtering...\n")
prof_all <- NULL
for (grid_id in 1:grids_n) {

  prof_cur_filter <- func_process_profile(vel_lbranch[,c(2,grid_id+2)])[,c(1,2)]
  prof_cur_filter$name <- paste0(substr(names(vel_lbranch)[grid_id+2], 1,4), "-", substr(names(vel_lbranch)[grid_id+2], 12,15))
  prof_all <- rbind(prof_all, prof_cur_filter)

}




#### Compute uncertainties, as NMAD of velocity magnitude on stable terrain  -----------------------
stable_ras <- rast("ref/stable_mask_buf200_step8.tif")
dir_scenes <- "2-scenes_all_cropped_clahe/"
scenes_nmad <- rep(NA_real_, length(grids))
for (grid_id in 1:length(grids)) {

  cat(grid_id, "/", length(grids), "\n")

  # Load original (CLAHE'd) images, to find out
  # where they were saturated or missing data.
  # For that, we use focal() with sd(), if sd is 0
  # over a 3x3 window then we assume we are in a saturated area.
  grid_name <- vel_maps[grid_id]
  name1 <- substr(grid_name, 1, 10)
  name2 <- substr(grid_name, 12, 21)

  fn1 <- list.files(dir_scenes, pattern = paste0(name1, ".*\\.tif$"))
  fn2 <- list.files(dir_scenes, pattern = paste0(name2, ".*\\.tif$"))

  r1 <- rast(file.path(dir_scenes, fn1))
  r1_sd <- focal(r1, w = 3, fun="sd", expand = TRUE)
  r1_sd_mask <- resample(subst(r1_sd, 0, NA) > 0, stable_ras, method = "near")

  r2 <- rast(file.path(dir_scenes, fn2))
  r2_sd <- focal(r2, w = 3, fun="sd", expand = TRUE)
  r2_sd_mask <- resample(subst(r2_sd, 0, NA) > 0, stable_ras, method = "near")


  grid_cur_masked <- stable_ras * grids[[grid_id]] * r1_sd_mask * r2_sd_mask
  scenes_nmad[grid_id] <- mad(values(grid_cur_masked, mat=FALSE, na.rm = T))
}
scenes_nmad_vel <- scenes_nmad * 365.25 / as.numeric(vel_maps_df$dt)


# Create data frame with also error interval.
# At this stage, we also introduce the uncertainty due to seasonality:
# we compute it by assuming a seasonal variability by up to N % (N = 50 - conservative, from Suslov1980)
# compared to annual mean, then the uncertainty (*relative* to the mean annual value) is given by:
# (N/100) * |L| / (365+L), where L is the difference between 365 and the actual interval of the considered velocity.
# We combine this uncertainty with the NMAD, using standard error propagation.

prof_all_err <- prof_all
prof_all_err$vel_lower <- NA_real_
prof_all_err$vel_upper <- NA_real_
names_u <- sort(unique(prof_all$name))

seasonality_factor <- 50 / 100 # assume (very conservative!) 50 % seasonal variability
for (name_id in 1:length(names_u)) {

  ids_cur <- which(prof_all_err$name == names_u[name_id])
  dt_cur <- as.integer(vel_maps_df$dt[name_id])
  seasonality_err_cur <- prof_all_err$vel[ids_cur] * seasonality_factor * abs(dt_cur - 365.25) / (dt_cur)
  
  
  prof_all_err$vel_lower[ids_cur] <- prof_all_err$vel[ids_cur] - sqrt(scenes_nmad_vel[name_id]^2 + seasonality_err_cur^2)
  prof_all_err$vel_upper[ids_cur] <- prof_all_err$vel[ids_cur] + sqrt(scenes_nmad_vel[name_id]^2 + seasonality_err_cur^2)

}



#### Velocity plot on left branch -----------------------------------------------------------------
palette_ext <- RColorBrewer::brewer.pal(11, "RdYlBu")[c(1:2, 4:5, 8:11)]

p_lines_lbranch <- ggplot(prof_all_err) +
  geom_ribbon(aes(x = dist_along/1000, ymin = vel_lower, ymax = vel_upper, fill = name), alpha = 0.15) +
  geom_line(aes(x = dist_along/1000, y = vel, color = name)) +#, linewidth = 1) +
  geom_line(data = subset(prof_all, name == '1996-1997'),
            aes(x = dist_along/1000, y = vel, color = name)) +#, linewidth = 1) + # Po plot 1996-1997 last, better readability.
  scale_color_manual(name = "Interval",
                     values = palette_ext) +
  scale_fill_manual(name = "Interval",
                     values = palette_ext) +
  scale_x_continuous(limits = c(profile_crop_upper, max(prof_all$dist_along))/1000,
                     expand = expansion(c(0.15,0.02),0)) +
  scale_y_continuous(limits = c(0,120),
                     expand = expansion(0,0),
                     breaks = seq(0,120,20),
                     oob = scales::oob_keep) +
  xlab("Distance along profile [km]") +
  ylab(bquote("Mean annual velocity [m yr"^"-1"*"]")) +
  geom_label(data = data.frame(x = c(profile_crop_upper, max(prof_all$dist_along)),
                               y = c(10,10),
                               label = c("2", "3")),
             aes(x = x/1000, y = y, label = label),
             fontface = "bold",
             size = 5,
             # size = unit(10, "cm"),
             label.padding = unit(0.25, "lines")) +
  theme_bw(base_size = 16) +
  theme(
        legend.position = c(0,0.8),
        legend.justification = c(0,1),
        legend.text = element_text(margin = margin(0,0,0.5,0,"lines")),
        legend.title = element_text(face = "bold"),
        plot.margin = margin(14,240,14,14,unit = "pt"))
# p_lines_lbranch
# ggsave("vel.pdf", plot = p_lines_lbranch, width = 12, height = 8)





#### Left branch: profile map ---------------------------------------------------------------------
north_arr <- readPNG("north2.png")
map_lbranch <- ggplot(prof_line_lbranch) +
  geom_spatraster(data = stretch(bg), maxcell = 1e7) +
  scale_fill_gradient2(low = "#000000", mid = "#DDDDDD", high = "#FFFFFF", midpoint = 127, guide = "none") +
  geom_sf(color = "#FF0000") +
  scale_x_continuous(limits = c(719350, 720600), expand = expansion(0,0)) +
  scale_y_continuous(limits = c(4389050, 4392550), expand = expansion(0,0)) +
  geom_label(data = data.frame(x=c(lbranch_line_coords[1,1],lbranch_line_coords[nrow(lbranch_line_coords),1]),
                               y=c(lbranch_line_coords[1,2],lbranch_line_coords[nrow(lbranch_line_coords),2]),
                               label = c("2","3")),
             aes(x = x, y = y, label = label),
             fontface = "bold",
             size = unit(10, "cm"),
             label.padding = unit(0.30, "lines")) +
  annotation_custom(grid::rasterGrob(north_arr,
                                     x = 0.992, y = 0.008, just = c('right', 'bottom'), width = unit(3.7, "cm"))) +
  ggsn::scalebar(x.min = 719350, x.max = 720600,
                 y.min = 4389050, y.max = 4392550,
                 st.dist = 0.02,
                 st.size = unit(4.2, "cm"),
                 location = "topright",
                 anchor = c(x = 720460, y = 4389290),
                 height = 0.01,
                 transform = FALSE,
                 dist = 250,
                 dist_unit = "m") +
  coord_sf(datum = 32642) +
  theme_void(base_size = 18) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0,"pt"),
        plot.margin = margin(0,0,0,0,unit = "pt"))
# map_lbranch
# ggsave("map.pdf", map_lbranch, width = 10*0.3571429, height = 10)

pl <- ggdraw(p_lines_lbranch) +
  draw_plot(map_lbranch, x = 0.425, y = 0.125, height = 0.85, hjust = 0, vjust = 0)
mult <- 0.44
ggsave(paste0("abramov_vel_w128_step8_mag_clahe", version_suffix, ".png"), plot = pl, width = 50*mult, height = 20*mult)
