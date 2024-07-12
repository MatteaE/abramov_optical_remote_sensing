# Functions for annual profile plots.
# First a function to plot velocities with lines,
# then a function to plot velocities as a ribbon (for both tracks together),
# then a function to plot the profile map.


# Function to plot annual velocity profiles.
func_plot_annual_profiles <- function(vel_melt,
                                      years,
                                      labels) {
  
  palette_ext <- RColorBrewer::brewer.pal(11, "RdYlBu")[c(1, 3:4, 9, 10:11)]
  
  pl <- ggplot(vel_melt) +
    geom_ribbon(aes(x = dist_along/1000,
                    ymin = value - vel_uncertainty, ymax = value + vel_uncertainty,
                    fill = variable),
                alpha = 0.15) +
    geom_line(aes(x = dist_along/1000, y = value, color = variable)) +
    scale_color_manual(name = "Interval",
                       values = palette_ext,
                       labels = paste0(as.character(years-1), "-", as.character(years))) +
    scale_fill_manual(name = "Interval",
                       values = palette_ext,
                       labels = paste0(as.character(years-1), "-", as.character(years))) +
    scale_linetype_discrete(name = "Orbital track") +
    scale_y_continuous(limits = c(0,60), expand = expansion(0,0)) +
    xlab("Distance along profile [km]") +
    ylab(bquote("Mean annual velocity [m yr"^"-1"*"]")) +
    geom_label(data = data.frame(x = c(0.0, max(vel_melt$dist_along)),
                                 y = c(4,4),
                                 label = labels),
               aes(x = x/1000, y = y, label = label),
               size = 5,
               fontface = "bold",
               label.padding = unit(0.25, "lines")) +
    theme_bw(base_size = 16)
  
  return(pl)
  
}


# Function to plot annual velocity profiles from both orbital tracks, as a ribbon.
func_plot_annual_profiles_ribbon <- function(vel_minmax_melt,
                                             years,
                                             labels) {
  
  
  
  
  pl <- ggplot(vel_minmax_melt) +
    geom_ribbon(aes(x = dist_along, ymin = value_min, ymax = value_max, fill = variable), alpha = 0.5) +
    scale_fill_discrete(name = "Year", labels = as.character(years)) +
    scale_linetype_discrete(name = "Orbital track") +
    scale_y_continuous(limits = c(0,70), expand = expansion(0,0)) +
    xlab("Distance along profile [m]") +
    ylab(bquote("Mean annual velocity [m yr"^"-1"*"]")) +
    geom_label(data = data.frame(x = c(0.0, max(vel_minmax_melt$dist_along)),
                                 y = c(4,4),
                                 label = labels),
               aes(x = x, y = y, label = label),
               size = 5,
               fontface = "bold") +
    theme_bw(base_size = 16) +
    theme(axis.title.y = element_markdown())
  
  return(pl)
  
}


# Function to plot single-profile maps.
func_plot_profile_map <- function(prof_line,
                                  prof_line_coords,
                                  extent,
                                  labels) {
  
  pl_map <- ggplot(prof_line) +
    geom_spatraster_rgb(data = bg, max_col_value = 255, maxcell = 1e6) +
    geom_sf(color = "#FF0000") +
    scale_x_continuous(limits = extent[1:2], expand = expansion(0,0)) +
    scale_y_continuous(limits = extent[3:4], expand = expansion(0,0)) +
    geom_label(data = data.frame(x=c(prof_line_coords[1,1],prof_line_coords[nrow(prof_line_coords),1]),
                                 y=c(prof_line_coords[1,2],prof_line_coords[nrow(prof_line_coords),2]),
                                 label = labels),
               aes(x = x, y = y, label = label),
               fontface = "bold",
               size = 2,
               label.padding = unit(0.15, "lines")) +
    annotation_scale(location = "br", height = unit(0.1, "cm")) +
    annotation_north_arrow(location = "br", style = north_arrow_minimal) +
    coord_sf(datum = 32642) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.length = unit(0,"pt"),
          plot.margin = margin(0,0,0,0))
  
  return(pl_map)
  
}
