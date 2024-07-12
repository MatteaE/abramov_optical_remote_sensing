# This code selects pictures to feed to the cosicorr-annual glacier velocity pipeline.
# It separates the two orbital tracks which cover Abramov glacier - thus,
# it should work with any site covered by two tracks, with minimal modification.
# Input: a folder with all scenes together
# Output: a directory tree with scenes properly selected and categorized.
# Output structure:
# 1-sel
#   path<A,B>
#     <2017-2023>
#       b<1,2,3,4>

# Overview:
# . Loop on orbital path (A and B)
# . . Loop on year (1 to 7, i.e. 2016/17 to 2022/23)
# . . . Loop on band (1,2,3,4, corresponding to R, G, B and NIR)


library(terra)
library(tools)

repodir        <- "/PATH/TO/ALL/Sentinel-2/SCENES"
targetdir_base <- paste0("PATH/TO/TREE/OF/SELECTED/SCENES/1-sel")

years_target <- 2017:2023
years_n <- length(years_target)

dates_incr <- 5 # Revisit interval along a same path [days].


orb_path_names <- c("pathA", "pathB")
orb_path_basedates <- as.Date(c("2020/09/05", "2020/09/08")) # These two dates are representative of the two orbital tracks.


dir.create(targetdir_base, recursive = TRUE)

for (orb_path_id in 1:2) {
  
  orb_path_cur <- orb_path_names[orb_path_id]
  basedate_cur <- orb_path_basedates[orb_path_id]
  
  cat("===== Orbital path", orb_path_cur, "=====\n")
  
  dir.create(file.path(targetdir_base, orb_path_cur))
  
  for (year_id in 1:years_n) {
    
    year_cur <- years_target[year_id]
    
    cat("  Year", year_cur, "\n")
    
    dir.create(file.path(targetdir_base, orb_path_cur, year_cur))
    
    sapply(file.path(targetdir_base, orb_path_cur, year_cur, paste0("b", 1:4)),
           dir.create)
    
    
    
    # Compute acceptable dates.
    date_min <- as.Date(paste0(year_cur - 1, "/07/01"))
    date_max <- as.Date(paste0(year_cur,     "/09/30"))
    date_start <- basedate_cur - floor((basedate_cur - date_min)/dates_incr)*dates_incr
    dates_sel <- seq.Date(date_start, date_max, dates_incr)
    
    filenames_sel <- paste0(dates_sel, ".tif")
    filenames_all <- list.files(repodir, pattern = "\\.tif$")
    
    filenames_match <- na.omit(match(filenames_sel, filenames_all))
    
    filenames_n <- length(filenames_match)
    if (filenames_n > 0) {
      
      for (file_id in 1:filenames_n) {
        
        filename_orig <- filenames_all[filenames_match[file_id]]
        full_ras      <- rast(file.path(repodir, filename_orig))
        
        # Bands 1-4 are R, G, B, NIR: all MSI bands at 10 m.
        for (band_id in 1:4) {
          
          
          targetdir <- file.path(targetdir_base, orb_path_cur, year_cur, paste0("b", band_id))
          filename_target <- format(as.Date(file_path_sans_ext(filename_orig)), paste0("%Y%m%d_b", band_id, ".tif"))
          writeRaster(full_ras[[band_id]],
                      file.path(targetdir, filename_target))
        } # End loop on bands
      } # End loop on found files of interest
    } # End if found more than 0 files
  } # End loop on years
} # End loop on orbital paths
