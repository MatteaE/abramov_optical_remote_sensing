# Here we define a stable terrain mask based on
# Copernicus surface slope and rasterized RGI6.0.

library(terra)

slope_lim <- 30 # Maximum slope to consider non-glacierized terrain stable [degrees].

rgi <- vect("/PATH/TO/RGI/VECTOR/FILE")
dem <- rast("/PATH/TO/DEM/GRID.tif")

rast_blueprint <- rast("/PATH/TO/ONE/COSICORR/CORRELATION/RESULT.tif")
dem_extr <- project(dem, rast_blueprint, method = "cubic")

sl <- terrain(dem_extr, "slope", unit = "degrees")

rgi_rast <- rasterize(rgi, dem_extr)

stable_mask <- (is.na(rgi_rast)) * (sl < slope_lim)

writeRaster(stable_mask, "ref/stable_mask.tif", overwrite = TRUE)
