# This script performs a rectilinear box method analysis of glacier terminus changes.
# It expects one vector file per date; each vector file can have the digitized terminus of any number of glaciers (one entity per glacier, distinguished by attribute "id").
# In the Abramov case we have used one folder per sensor type and one vector file per acquisition of that sensor type.

# Algorithm:
# Open all termini vector files
# For each glacier:
# Select its box
# Select all its outlines from the scene-wise termini gpkgs
# For each outline:
# Split the box in two polygons using the outline (lwgeom::st_split)
# Compute (zonal) mean DEM height over each of the two polygons
# Take the polygon with the higher mean height (it is the glacierized one)
# Divide its area by the length of the shorter of the two sides of the box rectangle
# Save this as its relative length

# We put everything in a list called gl_termini_all.
# This is a list of lists (one sub-list per glacier).
# Each sub-list has the glacier id and a data.frame with date and length.

setwd("./")

library(cowplot)
library(exactextractr) # for exact_extract from the DEM
library(grid)
library(ggplot2)
library(lwgeom)
library(sf)
library(terra)

dir.create("plots/by_glacier", recursive = TRUE)  # Here we put one plot per glacier id - only the time series of length
dir.create("plots/by_glacier_combined")           # Here we put one plot per glacier id - two panels: time series of length, and hypsometric distribution

# These are the paths to all folders containing vector files of terminus positions.
termini_dirs <- c("")


box_fp <- "/PATH/TO/VECTOR/FILE/OF/RECTILINEAR/BOXES" # This vector file must have an "id" field to match rectilinear boxes to terminus positions.
dem_fp <- "/PATH/TO/DEM/GEOTIFF" # DEM used to find the upper/glacierized and lower/unglacierized parts of each rectilinear box.

date_ref_for_plot <- "" # YYYY/MM/DD date of the zero-length reference, for the plots.

#### Load all vector files to a list --------------------------------------------------------------
# We use GPKGs, easy to adapt as needed.
termini_all <- list()
for (termini_dir_cur in termini_dirs) {
  lf <- list.files(termini_dir_cur, pattern = "\\.gpkg$")
  lfn <- length(lf)
  
  for (i in 1:lfn) {
    fn_cur <- lf[i]
    date_cur_str <- substr(fn_cur, 1, 10)
    # date_cur <- as.Date(date_cur_str)
    termini_date_cur <- st_read(file.path(termini_dir_cur, fn_cur),quiet = TRUE)
    termini_all[[date_cur_str]] <- termini_date_cur
  }
}

termini_dates_n <- length(termini_all)

box <- st_read(box_fp, quiet = T)
dem <- rast(dem_fp)


#### Define indices of glaciers -------------------------------------------------------------------
# We exclude two unmapped glaciers ending in rock glaciers,
# and we don't include the indices of the lines used for uncertainty (101 to 105).
glacier_ids_all     <- 1    # Values of the "id" field to be considered.
glacier_discard_ids <- c()  # Additional values to exclude

glacier_ids <- setdiff(glacier_ids_all, glacier_discard_ids)


#### Loop on the glaciers -------------------------------------------------------------------------
gl_termini_all <- list() # Here we will put one list per glacier, with glacier id and a df with date and length_rel.
for (gl_id_cur_id in 1:length(glacier_ids)) {
  
  gl_id_cur <- glacier_ids[gl_id_cur_id]
  
  box_cur <- box[which(box$id == gl_id_cur),]
  # Compute width of the box. It is the length of the shorter side of the rectangle.
  # To get it, we extract the polygon points (5 of them: it is closed) and compute their distance.
  box_cur_width <- min(st_distance(st_cast(box_cur, "POINT")[1:4,], st_cast(box_cur, "POINT")[2:5,], by_element = TRUE))
  
  gl_cur_dates_str <- character(0)
  gl_cur_lengths <- numeric(0)
  
  #### . Loop: find all termini for the current glacier -------------------------------------------
  for (terminus_date_id in 1:termini_dates_n) {
    termini_date_cur <- termini_all[[terminus_date_id]]
    date_cur_str <- names(termini_all)[terminus_date_id]
    
    gl_in_date_cur <- which(termini_date_cur$id == gl_id_cur)
    
    if (length(gl_in_date_cur) == 1) {
      cat("Found terminus information for glacier", gl_id_cur, "at date", date_cur_str, "\n")
      
      gl_cur_dates_str <- c(gl_cur_dates_str, date_cur_str)
      
      box_split <- st_collection_extract(st_split(box_cur, termini_date_cur[gl_in_date_cur,]))
      
      split_n <- nrow(box_split)
      if (split_n != 2) {
        stop("error splitting the box at date ", date_cur_str, " -- glacier id ", gl_id_cur)
      }
      
      box_mean_eles <- as.numeric(exact_extract(dem, box_split, fun = "mean"))
      
      box_gl_id <- which.max(box_mean_eles) # This is the higher of the two parts of the split box, therefore it contains the glacier.
      box_gl_area <- st_area(box_split[box_gl_id,])
      
      #### . . Apply box method here: get relative length of current terminus ---------------------
      gl_length <- box_gl_area / box_cur_width
      gl_cur_lengths <- c(gl_cur_lengths, gl_length)
      
      
      
    } else if (length(gl_in_date_cur) > 1) {
      stop("found more than one terminus corresponding to the current glacier. There is a mistake in the termini from ", date_cur_str, "\n")
    } else {
      cat("No terminus for glacier", gl_id_cur, "at date", date_cur_str, "\n")
    }
    
  } # End iterate on terminus collections.
  
  
  #### . Sort and process lengths to common reference ---------------------------------------------
  gl_cur_dates <- as.Date(gl_cur_dates_str)
  gl_df <- data.frame(date = gl_cur_dates,
                      length = gl_cur_lengths)
  
  # Sort by increasing date.
  gl_df <- gl_df[sort.int(gl_df$date, index.return = TRUE)$ix,]
  
  # Compute relative length change. For glaciers which do not have an outline
  # at the date of zero reference, we use the closest later date.
  id_reference <- which.min(abs((gl_df$date - as.Date(date_ref_for_plot)) / (gl_df$date >= as.Date(date_ref_for_plot))))
  gl_df$length_rel <- gl_df$length - gl_df$length[id_reference]
  
  
  gl_df$gl_id <- gl_id_cur
  gl_termini_all[[gl_id_cur_id]] <- list(gl_id = gl_id_cur,
                                         gl_df = gl_df)
  
  
} # End iterate on glacier ids.

gl_df_all <- do.call(rbind, sapply(gl_termini_all, "[", "gl_df"))

#### Check how many glaciers have a terminus available for each date ------------------------------
dates_uq <- unique(gl_df_all$date)
dates_uq_n <- length(dates_uq)
scene_count_per_date <- integer(dates_uq_n)
for (date_id in 1:dates_uq_n) {
  scene_count_per_date[date_id] <- length(which(gl_df_all$date == dates_uq[date_id]))
}



#### Plot: each glacier length individually -------------------------------------------------------
pl_lengths <- list()
for (gl_id_id in 1:length(glacier_ids)) {
  gl_id <- glacier_ids[gl_id_id]
  pl_lengths[[gl_id_id]] <- ggplot(gl_df_all[gl_df_all$gl_id == gl_id,]) +
    geom_line(aes(x = date, y = length_rel, color = as.factor(gl_id)), linewidth = 0.5) +
    geom_point(aes(x = date, y = length_rel), color = "black", size = 0.2) +
    annotation_custom(grobTree(textGrob(as.character(gl_id), x=-0.13,  y=-0.13, hjust=1, vjust = 0.5,
                                        gp=gpar(col="black", fontface="bold", fontsize = 16)))) +
    scale_x_date(breaks = seq.Date(as.Date("1940/01/01"), as.Date("2030/01/01"), "40 years"),
                 limits = c(as.Date("1968/01/01"), as.Date("2024/01/01")),
                 date_labels = "%Y",
                 expand = expansion(0,0)) +
    scale_y_continuous(n.breaks = 3) +
    coord_cartesian(clip = "off") +
    # xlab("Date") +
    ylab("Length [m]") +
    theme_bw(base_size = 16) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(1,10,7,1, "pt"))
  ggsave(paste0("plots/by_glacier/", gl_id, ".png"), width = 2, height = 2)
}



#### Plot: altitudinal distribution of each glacier -----------------------------------------------
# We put it on top of the altitudinal distribution of all glaciers together.
dem_32642 <- rast(dem_fp)
cellsize <- xres(dem_32642)*yres(dem_32642)/1e6
outl_all <- st_read("/PATH/TO/FULL/GLACIER/OUTLINES/VECTOR/FILE", quiet = T)
outl_sel <- outl_all[outl_all$id_abra_neigh %in% glacier_ids,]
eles <- exact_extract(dem_32642, outl_sel)
names(eles) <- outl_sel$id_abra_neigh

eles_all <- do.call(rbind, eles)

ele_range <- range(eles_all$value)
eleband_size <- 20
elebands_lower <- seq(floor(ele_range[1]/eleband_size)*eleband_size, floor(ele_range[2]/eleband_size)*eleband_size, eleband_size)
elebands_empty <- data.frame(lower = elebands_lower,
                             upper = elebands_lower + eleband_size,
                             midpoint = elebands_lower + eleband_size / 2,
                             area = NA_real_)
elebands <- elebands_empty
for (eleband_id in 1:nrow(elebands)) {
  cat(eleband_id, "\n")
  ncells_cur <- sum(eles_all$coverage_fraction[(eles_all$value >= elebands$lower[eleband_id]) &
                                                 (eles_all$value < elebands$upper[eleband_id])])
  elebands$area[eleband_id] <- ncells_cur * cellsize
  
}
elebands$area_rel <- elebands$area / sum(elebands$area)

# Now compute glacier-wise elevation distribution.
pl_eledists <- list()
for (gl_id_id in 1:length(glacier_ids)) {
  gl_id <- glacier_ids[gl_id_id]
  elebands_gl <- elebands_empty
  eles_cur <- eles[[as.character(gl_id)]]
  
  for (eleband_id in 1:nrow(elebands_gl)) {
    ncells_cur <- sum(eles_cur$coverage_fraction[(eles_cur$value >= elebands_gl$lower[eleband_id]) &
                                                   (eles_cur$value < elebands_gl$upper[eleband_id])])
    elebands_gl$area[eleband_id] <- ncells_cur * cellsize
  }
  elebands_gl$area_rel <- elebands_gl$area / sum(elebands_gl$area)
  pl_eledists[[gl_id_id]] <- ggplot(elebands) +
    geom_col(aes(x = midpoint, y = area_rel), width = eleband_size, fill = "#BBBBEE") +
    geom_col(data = elebands_gl, aes(x = midpoint, y = area_rel), width = eleband_size, fill = "#EA4444", alpha = 0.5) +
    coord_flip() +
    xlab("Elevation [m asl]") +
    ylab("Area fraction [-]") +
    scale_x_continuous(n.breaks = 3) +
    scale_y_continuous(n.breaks = 2) +
    theme_bw()
  # ggsave(paste0("plots/by_glacier/", gl_id, "_eledist.pdf"), width = 8, height = 5)
}



#### Plot: multi-panel with length change and elevation distribution ------------------------------
for (gl_id_id in 1:length(glacier_ids)) {
  gl_id <- glacier_ids[gl_id_id]
  plot_grid(pl_lengths[[gl_id_id]] + annotation_custom(grobTree(textGrob(as.character(gl_id), x=-0.23,  y=-0.17, hjust=1, vjust = 0.5,
                                                                         gp=gpar(col="black", fontface="bold", fontsize = 16)))),
            pl_eledists[[gl_id_id]], align= "hv", axis = "tlr", greedy = T, ncol = 2)
  ggsave(paste0("plots/by_glacier_combined/", gl_id, ".pdf"), width = 4, height = 2)
}

