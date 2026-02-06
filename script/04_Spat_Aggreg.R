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
dir$sf <- here(dir$data, "shapefiles")
dir$final <- here(dir$data, "final")

# Create non existing directories
lapply(dir, function(i) dir.create(i, recursive = T, showWarnings = F))

# Execute script "00_Functions.R" to load functions
source(here(dir$script, "00_Functions.R"))

#===============================================================================
# 1)  Aggregate at the commune level ------
#===============================================================================

# Function for spatial aggregation from "00_Functions.R"
spatial_agg_era5(
  source_dir = dir$prepared,
  start_year = "2006",
  end_year = "2020",
  start_reference_year = "1971",
  end_reference_year = "2000"
                 )


