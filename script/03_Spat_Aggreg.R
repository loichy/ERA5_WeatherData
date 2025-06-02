#===============================================================================
# Description: (De)-aggregate weather data at a communal level using France'  
# official communal shapefile
#===============================================================================

# Clean memory 
rm(list=ls())
gc()

# Load package
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ecmwfr, tidyverse, dplyr, terra, maps, here, ncdf4, raster, climate, devtools, lubridate,
               sf, sp, rnaturalearth, Matrix, readr, zoo, gstat)

# List directories 
dir <- list()
dir$root <- here()
dir$data <- here(dir$root, "data") # Folder for data
dir$source <- here(dir$data, "source") # Folder for original data files
dir$prepared <- here(dir$data, "prepared") # Folder for prepared/created data files
dir$script <- here(dir$root, "script") # Folder for code analyses
dir$output <- here(dir$root, "output") # Folder for created figures

# Create non existing directories
lapply(dir, function(i) dir.create(i, recursive = T, showWarnings = F))

# Execute script "00_Functions.R" to load functions
source(here(dir$script, "00_Functions.R"))

#===============================================================================
# 1)  Aggregate using commune shapefile ------
#===============================================================================

# Function for spatial aggregation
spatial_agg_era5(source_dir = dir$prepared)

# Function for correcting WSD and CSD columns to integer values after interpolation
correct_interpol <- function(source_dir) {
  files <- list.files(source_dir, full.names = TRUE, pattern = "\\.rds$")
  files <- files[basename(files) != "normals.rds"]
  
  for (file in files) {
    final_centroid_interp <- readRDS(file)
    cols_to_round <- c("warm_spell", "warm_spell_index", "cold_spell", "cold_spell_index")
    
    final_centroid_interp[cols_to_round] <- lapply(final_centroid_interp[cols_to_round], function(col) {
      if (is.numeric(col)) {
        return(round(col))
      } else {
        return(col)
      }
    })
    saveRDS(final_centroid_interp, file) 
  }
}

correct_interpol(source_dir = dir$prepared)



