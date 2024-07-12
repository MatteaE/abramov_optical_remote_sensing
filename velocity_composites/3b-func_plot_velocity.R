# Function to plot velocity (color-coded + arrows).
func_plot_velocity <- function(vel_ew,
                               vel_ns,
                               vel_mag) {
  
  
  #### . Prepare data for velocity plots ------------------------------------------------------------
  # To plot arrows and/or streamlines.
  vel_crds <- crds(vel_ew, na.rm = FALSE)
  vel_df_streamlines <- data.frame(x = vel_crds[,1],
                                   y = vel_crds[,2],
                                   dx = values(vel_ew, mat = FALSE),
                                   dy = values(vel_ns, mat = FALSE),
                                   vel = values(vel_mag, mat = FALSE))
  
  
  
  
  # after_stat() to retrieve the layer without using its name
  # (it is "median" for band plots "sum" for composite).
  pl_vel <- ggplot() +
    geom_spatraster(data = vel_mag, aes(fill = after_stat(value))) +
    metR::geom_arrow(data = vel_df_streamlines,             # Arrow heads do also scale with velocity.
                     aes(x = x, y = y, dx = dx, dy = dy),
                     skip.x = 3, skip.y = 3,
                     arrow.angle = 30, arrow.type = "open", arrow.length = unit(4, "pt"),
                     pivot = 0, preserve.dir = TRUE, direction = "ccw", size= .6) +
    scale_fill_distiller(name = "Median annual\nvelocity [m/yr]",
                         type = "seq", palette = "Reds", direction = 1,
                         oob = scales::squish,
                         limits = c(0,60)) +
    scale_mag(max = 80,
              max_size = 0.53,
              guide = "none") +
    coord_sf(datum = 32642,
             xlim = c(714100, 723860),
             ylim = c(4385100, 4393860),
             expand = FALSE) +
    xlab("UTM 42N Easting [m]") +
    ylab("UTM 42N Northing [m]") +
    theme_bw()
  
  
  return(pl_vel)
  
}
