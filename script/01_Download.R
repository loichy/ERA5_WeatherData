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

key_loic <- "bcf2827f-fb1f-4647-ac97-e52ba81be39c"  

# Key Inès: 73737aee-7063-4fb9-9548-2a983ffdbfbe
# Key Loic: bcf2827f-fb1f-4647-ac97-e52ba81be39c

#===============================================================================
# 2) Download weather data ------
#===============================================================================

# Define combinations of variables and statistics
params <- tribble(
  ~variable,        ~statistic,
  "2m_temperature", "daily_mean",
  "2m_temperature", "daily_minimum",
  "2m_temperature", "daily_maximum",
  "total_precipitation", "daily_sum"
)

# Define other shared parameters
start_date <- "2018-06-01"
end_date <- "2019-12-31"
output_dir <- dir$source

# Run download function for each combination
params |>
  filter(variable == "total_precipitation") %>% 
  purrr::pwalk(~ download_era5_ecmwfr(
    variable = ..1,
    statistic = ..2,
    start_date = start_date,
    end_date = end_date,
    output_dir = output_dir,
    key = key_loic
  ))

#===============================================================================
# 3) Download Land-Sea Mask ------
#===============================================================================

# Go to Climate Data Store (CDS) to download the Land-Sea Mask manually (not working with ecmwfr package)



