# Function to process a longitudinal velocity profile.
# . Prune points (so that there is only one point for each velocity retrieval, associated with the middle dist_along)
# . Remove NAs
# . Filter on longitudinal deformation rates (compression/extension)
# . Re-introduce NAs, for the plot: each profile should not have gaps within its
# dist_along larger than the diagonal of the cell of the velocity raster,
#   so we insert NAs in them.

# df_profile has dist_along and the velocity.
func_process_profile <- function(df_profile_raw) {
  
  df_profile_v2 <- df_profile_raw
  names(df_profile_v2) <- c("dist_along", "vel")
  
  vel_rle <- rle(df_profile_v2$vel)
  
  ids_rle_multiple <- which(vel_rle$lengths > 1)
  ids_rle_multiple_n <- length(ids_rle_multiple)
  
  if (ids_rle_multiple_n > 0) {
    
    for (m_id in 1:ids_rle_multiple_n) {
      
      m_id_cur <- ids_rle_multiple[m_id]
      m_id_start_id <- sum(vel_rle$lengths[1:(m_id_cur-1)]) + 1
      m_id_end_id <- m_id_start_id + vel_rle$lengths[m_id_cur] - 1
      
      vel_cur <- df_profile_v2$vel[m_id_start_id]
      # We pick the new representative dist_along as mean of
      # the first and last (to deal with non-equidistant points).
      dist_along_cur_mean <- mean(df_profile_v2$dist_along[c(m_id_start_id,m_id_end_id)])
      
      # Apply representative dist_along to entries with a same velocity.
      # Then we will remove duplicated lines in the data.frame.
      df_profile_v2$dist_along[m_id_start_id:m_id_end_id] <- dist_along_cur_mean
      
    }
  }
  
  #### Remove duplicate velocity retrievals -------------------------------------------------------
  # Like this we only have NAs and 1 velocity retrieval per velocity.
  ids_dup <- which(duplicated.data.frame(df_profile_v2))
  if (length(ids_dup) > 0) {
    df_profile_v3 <- df_profile_v2[-ids_dup,]
  } else {
    df_profile_v3 <- df_profile_v2
  }
  
  # Remove also NAs.
  ids_na <- which(is.na(df_profile_v3$vel))
  if (length(ids_na) > 0) {
    df_profile_v4 <- df_profile_v3[-ids_na,]
  } else {
    df_profile_v4 <- df_profile_v3
  }
  
  
  #### Now compute longitudinal strain rates, for filtering ---------------------------------------
  df_profile_v5 <- df_profile_v4
  df_profile_v5$dvel <- c(0.0, diff(df_profile_v5$vel))
  df_profile_v5$ddist <- c(0.0, diff(df_profile_v5$dist_along))
  df_profile_v5$strain_rate <- df_profile_v5$dvel / df_profile_v5$ddist
  df_profile_v5$strain_rate[1] <- 0.0 # It was 0.0 / 0.0 = NaN.
  
  strain_rate_threshold <- 0.1
  
  exit_logi <- all(abs(df_profile_v5$strain_rate) < strain_rate_threshold)
  while (exit_logi == FALSE) {
    
    # This is the index of the largest strain rate, which is measured between this point and the PREVIOUS one.
    id_strain_worst <- which.max(abs(df_profile_v5$strain_rate))
    
    # cat("Max strain rate:", round(df_profile_v5$strain_rate[id_strain_worst], 3), "\n")
    
    # Remove the point causing excess strain rate.
    id_to_remove_id <- which.min(df_profile_v5$vel[id_strain_worst + c(-1,0)])
    df_profile_v5 <- df_profile_v5[-(id_strain_worst + (id_to_remove_id - 2)),]
    
    # Recompute strain rates.
    df_profile_v5$dvel <- c(0.0, diff(df_profile_v5$vel))
    df_profile_v5$ddist <- c(0.0, diff(df_profile_v5$dist_along))
    df_profile_v5$strain_rate <- df_profile_v5$dvel / df_profile_v5$ddist
    df_profile_v5$strain_rate[1] <- 0.0 # It was 0.0 / 0.0 = NaN.
    
    exit_logi <- all(abs(df_profile_v5$strain_rate) < strain_rate_threshold)
    
  }
  
  df_profile_v6 <- df_profile_v5[,c(1,2,4)]
  
  
  # Find isolated points (farther from their neighbor than a threshold, on both sides)
  # and extend them to the cell size, to have them more visible on the plot.
  # We limit this extension to avoid overlaps (the min() and max() below),
  # and also to stay at or above dist_along = 0.0.
  ddist_thresh_detect <- 40   # Detect isolated points (farther than cell size on both sides)
  ddist_extend        <- 20   # Extend on both sides by maximum half the cell size.
  isolated_ids <- which((c(ddist_thresh_detect+1, df_profile_v6$ddist[2:nrow(df_profile_v6)]) > ddist_thresh_detect) &
                                (c(df_profile_v6$ddist[2:nrow(df_profile_v6)], ddist_thresh_detect+1) > ddist_thresh_detect))
  isolated_ids_n <- length(isolated_ids)
  
  if (isolated_ids_n > 0) {
    
    # We have to make a quirk in case the first or last points are isolated,
    # we can't just split the vectors and insert in their case.
    if (1 %in% isolated_ids) {
      
      vels <- c(rep(df_profile_v6$vel[1], 2), df_profile_v6$vel)
      dists_along <- c(max(0.0, df_profile_v6$dist_along[1] - ddist_extend),
                       df_profile_v6$dist_along[1],
                       min(df_profile_v6$dist_along[1] + ddist_extend, df_profile_v6$dist_along[2]),
                       df_profile_v6$dist_along[2:nrow(df_profile_v6)])
      ddists <- c(0.0, diff(dists_along))
      df_profile_v6 <- data.frame(dist_along = dists_along,
                                  vel = vels,
                                  ddist = ddists)
      isolated_ids <- which((c(ddist_thresh_detect+1, df_profile_v6$ddist[2:nrow(df_profile_v6)]) > ddist_thresh_detect) &
                              (c(df_profile_v6$ddist[2:nrow(df_profile_v6)], ddist_thresh_detect+1) > ddist_thresh_detect))
    }
    dfnr <- nrow(df_profile_v6)
    if (dfnr %in% isolated_ids) {
      
      vels <- c(df_profile_v6$vel, rep(df_profile_v6$vel[dfnr], 2))
      dists_along <- c(df_profile_v6$dist_along[1:(dfnr-1)],
                       max(df_profile_v6$dist_along[dfnr-1], df_profile_v6$dist_along[dfnr] - ddist_extend),
                       df_profile_v6$dist_along[dfnr],
                       df_profile_v6$dist_along[dfnr] + ddist_extend)
      ddists <- c(0.0, diff(dists_along))
      df_profile_v6 <- data.frame(dist_along = dists_along,
                                  vel = vels,
                                  ddist = ddists)
      isolated_ids <- which((c(ddist_thresh_detect+1, df_profile_v6$ddist[2:nrow(df_profile_v6)]) > ddist_thresh_detect) &
                              (c(df_profile_v6$ddist[2:nrow(df_profile_v6)], ddist_thresh_detect+1) > ddist_thresh_detect))
    }
    

    # Recompute number of isolated ids, in case something was changed at the edges.    
    isolated_ids_n <- length(isolated_ids)
    
    # Now loop to extend isolated points within each profile (not at the edges).
    while (isolated_ids_n > 0) {
      
      vels <- c(df_profile_v6$vel[1:(isolated_ids[1]-1)],
                rep(df_profile_v6$vel[isolated_ids[1]], 3),
                df_profile_v6$vel[(isolated_ids[1]+1):nrow(df_profile_v6)])
      dists_along <- c(df_profile_v6$dist_along[1:(isolated_ids[1]-1)],
                       max(df_profile_v6$dist_along[isolated_ids[1]-1], df_profile_v6$dist_along[isolated_ids[1]] - ddist_extend),
                       df_profile_v6$dist_along[isolated_ids[1]],
                       min(df_profile_v6$dist_along[isolated_ids[1]] + ddist_extend, df_profile_v6$dist_along[isolated_ids[1]+1]),
                       df_profile_v6$dist_along[(isolated_ids[1]+1):nrow(df_profile_v6)])
      ddists <- c(0.0, diff(dists_along))
      df_profile_v6 <- data.frame(dist_along = dists_along,
                                  vel = vels,
                                  ddist = ddists)
      isolated_ids <- which((c(ddist_thresh_detect+1, df_profile_v6$ddist[2:nrow(df_profile_v6)]) > ddist_thresh_detect) &
                              (c(df_profile_v6$ddist[2:nrow(df_profile_v6)], ddist_thresh_detect+1) > ddist_thresh_detect))
      isolated_ids_n <- length(isolated_ids)
      
    }
    
  }
  
  

  
  #### Now reinsert NAs, for the plot -------------------------------------------------------------
  # Longest allowed distance between two points to consider them consecutive.
  # All gaps with a distance greater than this one will get an NA velocity in the middle.
  ddist_thresh <- 40 * sqrt(2) 
  ddist_too_long_ids <- which(df_profile_v6$ddist > ddist_thresh)
  
  while (length(ddist_too_long_ids) > 0) {
    
    dists_along <- c(df_profile_v6$dist_along[1:(ddist_too_long_ids[1] - 1)],
                     NA_real_,
                     df_profile_v6$dist_along[ddist_too_long_ids[1]:nrow(df_profile_v6)])
    dists_along <- zoo::na.approx(dists_along) # The dist_along of the new point is linearly interpolated.
    
    # The velocity of the new point remains NA, it is used to break the line plot.
    vels <- c(df_profile_v6$vel[1:(ddist_too_long_ids[1] - 1)],
              NA_real_,
              df_profile_v6$vel[ddist_too_long_ids[1]:nrow(df_profile_v6)])
    
    ddists <- c(0.0, diff(dists_along))
    df_profile_v6 <- data.frame(dist_along = dists_along,
                                vel = vels,
                                ddist = ddists)
    ddist_too_long_ids <- which(df_profile_v6$ddist > ddist_thresh)
    
  }
  
  return(df_profile_v6)
  
}
