# This script processes the ASP tie points that we will
# feed into Cosi-Corr:

# - Convert the tie points from pixel coordinates to lat/lon
# - Check which tie points are too close to the border of either the ref or the raw image, reject them
# - Check which tie points fall inside a mask (of unstable ground), reject them
# - If there are more than <max_gcps> points, reject the excess ones (else CosiCorr becomes super slow refining each one). To do this: select the first point randomly, then the second is the one farthest away, the third maximizes average distance from the first two, and so on. This ensures the best coverage on the image, better than the old version of selecting points with a (reproducible) random sampling.
# - Plot all GCPs against image background, colored according to accepted/rejected(/excess)
# - Write final file with only accepted tie points

func_process_matches <- function(
    
  ref_p,      # Path to reference orthoimage (GeoTiff)
  raw_p,      # Path to raw image to be orthorectified (GeoTiff)
  tps_p,      # Path to text file with the tie points in image coordinates (x1 x2 y1 y2)
  mask_p,     # Path to vector file with the unstable terrain mask
  gsd,        # Ground sampling distance (20 or 10 or 5 m). Used to reject TPs based on proximity to image border, within the size of the correlation window (NOT the half-size, as seen in CosiCorr's patches plots under RSM_Refinement/).
  max_gcps    # Keep only this many GCPs.

  ) {
  
  
  corr_win_mult <- 2 # We take this as safety factor for the margin at the borders of the image.
  if (gsd == 20) {
    corr_win_margin <- 32 * corr_win_mult
  } else {
    corr_win_margin <- 64 * corr_win_mult
  }
  
  tps   <- read.table(tps_p, header = F)
  names(tps) <- c("x1", "y1", "x2", "y2")
  
  
  ref_r <- rast(ref_p)
  
  # As of Feb 2024, terra complains that raw_r is a rotated raster, and refuses to load data from it.
  # Since here we are not interested in the georeferencing information of raw_r (only the ncol and nrow,
  # and we plot raw_r in pixel coordinates), we make a copy of raw_r where we discard all georeference.
  raw_p_nogeoref <- paste0(tools::file_path_sans_ext(raw_p), "_nogeoref.TIF")
  file.copy(raw_p, raw_p_nogeoref)
  system2("gdal_edit.py",
          args = c("-unsetgt", raw_p_nogeoref))
  raw_r <- rast(raw_p_nogeoref, lyrs = 1) # Only first layer in case there are several.
  crs(raw_r) <- ""
  ext(raw_r) <- ext(0, ncol(raw_r), 0, nrow(raw_r))
  
  
  # We have verified that these below correspond exactly to the map coordinates computed by Cosi-Corr.
  gcps <- cbind(xFromCol(ref_r, 1) + (tps$x1)*xres(ref_r),
                yFromRow(ref_r, 1) - (tps$y1)*yres(ref_r))
  gcps_vec <- vect(data.frame(x  = gcps[,1],
                              y  = gcps[,2],
                              id = 1:nrow(gcps)),
                   geom = c("x", "y"),
                   crs = "EPSG:32642")
  
  #### GCPs selection 1: which GCPs are on the border of the ref or raw image? ----------------------
  # Compute a SpatVector with 1 square polygon per GCP,
  # same size as the correlation window size. We use this to
  # check whether cosi-corr refinement is at risk of including
  # NA values (of the reference image) within the correlation.
  # We also add 1 pixel margin on all sides of the square
  # since the squares are not aligned with the grids.
  gcps_xmin <- gcps[,1] - ((corr_win_margin+1) * gsd)
  gcps_xmax <- gcps[,1] + ((corr_win_margin+1) * gsd)
  gcps_ymin <- gcps[,2] - ((corr_win_margin+1) * gsd)
  gcps_ymax <- gcps[,2] + ((corr_win_margin+1) * gsd)
  gcps_squares_coords <- cbind(gcps_xmin, gcps_xmax, gcps_ymin, gcps_ymax)
  gcps_n <- nrow(gcps_squares_coords)
  
  sq_list <- list()
  for (gcp_id in 1:gcps_n) {
    
    sq_list[[gcp_id]] <- rbind(c(gcps_xmin[gcp_id], gcps_ymin[gcp_id]),
                               c(gcps_xmax[gcp_id], gcps_ymin[gcp_id]),
                               c(gcps_xmax[gcp_id], gcps_ymax[gcp_id]),
                               c(gcps_xmin[gcp_id], gcps_ymax[gcp_id]))
  }
  gcps_squares <- vect(sq_list, "polygons")
  crs(gcps_squares) <- "EPSG:32642"
  
  # This computes the sum of the raster within each square.
  # If the square goes outside the data extent, the sum will be NaN.
  square_sums <- terra::zonal(ref_r, gcps_squares, fun = sum)
  
  # These are the GCPs whose correlation window includes nodata regions on the orthorectified reference.
  square_bad_ids_ref <- which(is.na(square_sums[,1]))
  
  # Now we look for GCPs whose correlation window includes nodata regions on the raw image.
  # These are easier - we simply have to look at squares in pixel (TPs) coordinates.
  # But we have to add some more margin since CosiCorr appears to be discarding a collar
  # around the raw image, for some reason!! So we double again the margin here.
  square_bad_ids_raw <- which(((tps$x2 + corr_win_margin*2 + 1) > ncol(raw_r)) |
                                ((tps$x2 - corr_win_margin*2 - 1) < 1) |
                                ((tps$y2 + corr_win_margin*2 + 1) > nrow(raw_r)) |
                                ((tps$y2 - corr_win_margin*2 - 1) < 1))
  
  gcps_border_bad_ids <- unique(c(square_bad_ids_ref, square_bad_ids_raw))
  
  # gcps_squares_fitered <- gcps_squares[setdiff(1:gcps_n, gcps_border_bad_ids),]
  # writeVector(gcps_squares_fitered, "gcps_filtered.gpkg")
  
  
  #### GCPs selection 2: which GCPs are on unstable terrain? ----------------------------------------
  unstable_mask <- vect(mask_p)
  gcps_vec_stable <- mask(gcps_vec,
                          unstable_mask,
                          inverse = TRUE)
  gcps_unstable_bad_ids <- setdiff(1:nrow(gcps), gcps_vec_stable$id)
  
  
  # The "ok" GCPs: on stable terrain, and not close to the scenes borders.
  gcps_vec_ok <- gcps_vec_stable[which(!(gcps_vec_stable$id %in% gcps_border_bad_ids)),]
  
  
  gcps_all_bad_ids <- c(gcps_border_bad_ids, gcps_unstable_bad_ids)
  gcps_vec_notok <- gcps_vec[gcps_all_bad_ids,]
  
  # For second plot, in raw image pixel coordinates.
  df_raw <- data.frame(x = tps[,3], y = nrow(raw_r) - tps[,4] + 1)
  df_raw_ok <- df_raw[gcps_vec_ok$id,]
  df_raw_notok <- df_raw[setdiff(1:nrow(gcps), gcps_vec_ok$id),]
  
  gcps_ok_n <- nrow(gcps_vec_ok)
  
  #### Remove excess GCPs ---------------------------------------------------------------------------
  # Select the first point randomly, then the second
  # is the one farthest away, the third maximizes average
  # distance from the first two, and so on.
  # Like this we have to compute the distance matrix
  # just once.
  if (gcps_ok_n > max_gcps) {
    
    dist_all <- terra::distance(gcps_vec_ok)
    dist_all_mat <- as(dist_all, "matrix")
    
    set.seed(1)
    ids_keep <- sample(1:gcps_ok_n, 1, replace = FALSE)
    
    # We also take the 2nd gcp manually before the loop,
    # to avoid selecting a single column from the distance
    # matrix later (else there is a type error).
    id_p2 <- which.max(dist_all_mat[,ids_keep])
    ids_keep <- c(ids_keep, id_p2)
    
    
    for (id_cur in 3:max_gcps) {
      
      # For each gcp, compute average distance to the set of points selected so far.
      # We will select the gcp for which this distance is highest.
      # NOTE: in some rare cases, it might be possible that the selected GCP was
      # already selected. So, we also compute an extra column vector (vec_gcp_available_logi)
      # which is 0 for all GCPs which have already being selected, and 1 for all
      # GCPs which are still available.
      # Then we can multiply it by the average distance vector to be sure that
      # we take a still unused GCP.
      mat_ids_keep <- dist_all_mat[,ids_keep]
      mat_ids_keep_rowMins <- rowMins(mat_ids_keep)
      vec_gcp_available_logi <- rep(1, gcps_ok_n)
      vec_gcp_available_logi[ids_keep] <- 0
      mat_ids_keep_rowMins_available <- mat_ids_keep_rowMins * vec_gcp_available_logi
      
      id_next <- which.max(mat_ids_keep_rowMins_available)
      
      ids_keep <- c(ids_keep, id_next)
    }
    
    ids_excess <- setdiff(1:gcps_ok_n, ids_keep)
    
    gcps_vec_excess <- gcps_vec_ok[ids_excess,]
    df_raw_excess <- df_raw_ok[ids_excess,]
    
    gcps_vec_ok <- gcps_vec_ok[ids_keep,]
    df_raw_ok <- df_raw_ok[ids_keep,]
  }
  
  
  tps_out <- tps[gcps_vec_ok$id,]
  
  # Plot reference orthoimage with selected and discarded tie points.
  pl_ref <- ggplot() +
    geom_spatraster(data = stretch(ref_r, histeq = TRUE),
                    maxcell = 2e6) +
    scale_fill_gradient(low = "#000000", high = "#FFFFFF", na.value = "#00000000",
                        oob = scales::squish,
                        guide = "none") +
    geom_spatvector(data = gcps_vec_notok, color = "red", shape = 4, size = 7, stroke = 1.5) +
    {if (gcps_ok_n > max_gcps) geom_spatvector(data = gcps_vec_excess, color = "yellow", shape = 4, size = 7, stroke = 1.5)} +
    geom_spatvector(data = gcps_vec_ok, color = "green", shape = 3, size = 10, stroke = 3) +
    coord_sf(datum = "EPSG:32642", expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0,0,0,0))
  
  
  # Same plot but for the raw image (using pixel coordinates).
  # margin_px <- 100 # Used to zoom in on the GCPs area only.
  pl_raw <- ggplot() +
    geom_spatraster(data = raw_r,
                    maxcell = 2e6) +
    scale_fill_gradient(low = "#000000", high = "#FFFFFF", na.value = "#00000000",
                        limits = c(-31,249),
                        oob = scales::squish,
                        guide = "none") +
    geom_point(data = df_raw_notok, aes(x = x, y = y), color = "red", shape = 4, size = 7, stroke = 1.5) +
    {if (gcps_ok_n > max_gcps) geom_point(data = df_raw_excess, aes(x = x, y = y), color = "yellow", shape = 4, size = 7, stroke = 1.5)} +
    geom_point(data = df_raw_ok, aes(x = x, y = y), color = "green", shape = 3, size = 10, stroke = 3) +
    coord_equal(expand = FALSE) +
    scale_x_continuous(#limits = c(min(df_raw$x - margin_px), max(df_raw$x + margin_px)), # Zoom in on the GCPs area only. 
                       expand = expansion(0,0)) +
    scale_y_continuous(#limits = c(min(df_raw$y - margin_px), max(df_raw$y + margin_px)), # Zoom in on the GCPs area only. 
                       expand = expansion(0,0)) +
    theme_void() +
    theme(plot.margin = margin(0,0,0,0))
  
  
  # Cleanup un-georeferenced raw_r that we used for plot.
  file.remove(raw_p_nogeoref)
  
  return(list(tps_out = tps_out,
              pl_ref = pl_ref,
              pl_raw = pl_raw
  ))
  
}
