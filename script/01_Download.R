#===============================================================================
# Description: Download Copernicus weather data using ecmwfr package
#===============================================================================

# Clean memory 
rm(list=ls())
gc()

# Load package
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ecmwfr, tidyverse, dplyr, terra, maps, here, ncdf4, raster, climate, devtools, lubridate,
               sf, sp, rnaturalearth, Matrix)

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
# 1) Set API key ------
#===============================================================================

# 1. Register to ECMWF to get an account and be able to request the API: https://www.ecmwf.int/

# 2. Get your API keys: https://cds.climate.copernicus.eu/how-to-api

# 3. Download weather data of interest

download_era5_ecmwfr(
  variable = "2m_temperature",
  start_date = "2023-01-01",
  end_date = "2023-02-28",
  statistic = "daily_mean",
  key = "73737aee-7063-4fb9-9548-2a983ffdbfbe",
  output_dir = dir$source
)




