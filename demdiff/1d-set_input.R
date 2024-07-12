# This file sets the input specs (paths and todo list) for the DEMdiff analysis.



#### Define dem_pool, in chronological order --------------------------------------------------------------------------
# These are simply the paths to all DEMs to be considered.
dem_pool <- normalizePath(c(""))


# Define dates of the DEMs in the above pool.
dem_dates <- as.Date(c(""))


#### Define outl_pool, in the order used for processing ---------------------------------------------------------------
# These are the path to the vector files of outlines used as aggregation polygons.
outl_pool <- normalizePath(c(""))



#### Define pool of grids of slope, aspect and curvature (max and min) ------------------------------------------------
# These are used for the analysis of terrain-dependent biases.
# We can have several in case we don't always work at the
# same resolution (e.g. keeping Pléiades at high resolution).
# All are selected using a same id in the todo_df: biascorr_topo_grids_id.
slope_pool  <- normalizePath(c(""))
aspect_pool <- normalizePath(c(""))
maxc_pool   <- normalizePath(c(""))
minc_pool   <- normalizePath(c(""))
bias_plots_ylim_multi_pool <- c() # Vertical limits of the plots in the bias analysis. One value per item of the previous pools.






#### Define pool of polygon masks of unstable terrain for coregistration ----------------------------------------------
# These are again paths to vector files.
# More than one in case unstable terrain changed a lot during the considered periods
# (e.g. use RGI for recent comparisons, but 1970s extents for old ones).
# NOTE: these outlines will be combined with the slopes and NMAD mask to compute the final unstable terrain.
unstable_outlines_pool <- normalizePath(c(""))




#### Declare and fill todo_df -----------------------------------------------------------------------------------------
todo_df <- data.frame(
  tba_dem_id            = integer(0),     # Id of the tba/second DEM within the DEM pool.
  ref_dem_id            = integer(0),     # Id of the reference/first DEM within the DEM pool.
  tba_dem_preproc       = integer(0),     # Flag for preprocessing of the tba/second DEM: 1 = crop/extend to extent, 2 = resample to blueprint, 3 = reproject to blueprint.
  ref_dem_preproc       = integer(0),     # Flag for preprocessing of the reference/first DEM: 1 = crop/extend to extent, 2 = resample to blueprint, 3 = reproject to blueprint.
  tba_coreg_flag        = integer(0),     # 0 = never coregister TBA, 1 = coregister TBA once (after preprocessing), 2 = coregister TBA twice (after preprocessing and again after bias correction; only makes sense if we do apply some bias correction, i.e. if any item in biascorr_flag is not 0).
  coreg_outl_id         = integer(0),     # Id of the unstable outlines file to load, for coreg. NOTE: these get combined with slopes and MAD mask.
  biascorr_flag         = character(0),   # Which bias corrections should be applied (vs simply quantified)? Binary string, e.g. "1110000". See 1e-biascorr_func.R.
  biascorr_topogrids_id = integer(0),     # Index of the slope/aspect/curvature grids from the respective pools. Also used for the ylim scale of the bias plots.
  biascorr_along_track  = character(0),   # This can be "SPOT", "ASTER", "Pléiades", or the character of a numeric. If it is the name of a satellite, it is used to compute the along-track angle (Pléiades has 0.0 since its undulations are North-South). If it is convertible to a numeric, it is assumed to already be the correct along-track angle (used e.g. for KH-9 which we don't compute explicitly). This angle is used to estimate/correct bias.
  outl_file_id          = integer(0),     # Index of the vector file with the polygon outlines to consider within this comparison. Used to select from outl_pool.
  outl_ids_filter       = character(0),   # Dash-separated indices of polygons to use for hypsometric filtering of the demdiff. As character! It can also be "all", then it uses every outline in the vector file.
  outl_ids_aggr         = character(0),   # Dash-separated indices of polygons over which we compute aggregated elevation change. As character! It can also be "all", then it uses every outline in the vector file.
  plot_scale_mult       = numeric(0)      # Multiplier of the color scale in plots of the DEM difference.
)



# Example of definition of a DEM comparison to be analyzed.
todo_df[ 1,] <- list(11,  8, 1, 1, 1, 2, "0000000", 1,     "SPOT",  1, "all", "1-8-400", 4)
