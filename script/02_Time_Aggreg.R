#===============================================================================
# Description: Delimit the raster borders to France and transform all files 
# into a data frame format to create time-aggregated indicators
#===============================================================================
ok
# Clean memory 
rm(list=ls())
gc()

# Load package
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ecmwfr, tidyverse, dplyr, purrr, terra, maps, here, ncdf4, raster, climate, devtools, lubridate,
               sf, sp, rnaturalearth, rnaturalearthdata, Matrix)

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
# 1) Prepare weather data into data frame format ------
#===============================================================================

df_list <- process_era5_files(source_dir = "data/source/")

#===============================================================================
# 2) Create time indices from monthly weather data ------
#===============================================================================

# For each month data frame in the list, create mean, min and max of month columns
df_list_time_agg <- lapply(df_list, function(df) weather_df_indice(df, prefixe_variable = "t2m", output_dir = dir$prepared))

#===============================================================================
# 3) Prepare data frame containing season normals that will be needed for WSDI
# and CSDI computations ------
#===============================================================================

# 1). Load raster of daily mean temperatures over 20 years in France
rast_normals <- readRDS(file = "data/source/MeanTemps_France_Daily_19812000.rds")

# 2). Clean raster data and transform into data frame format

# Get France borders with sf package
france <- ne_countries(scale = "medium", country = "France", returnclass = "sf")

# Harmonize projection
france_proj <- st_transform(france, crs(rast_normals))

# Crop data using France borders
cropped <- crop(rast_normals, vect(france_proj))
masked <- mask(cropped, vect(france_proj))

# Transform into data frame format
df_normals <- as.data.frame(masked, na.rm = TRUE)

# Extract xy coordinates to add to the final normals dataset
valid_pixel <- which(!is.na(values(masked)[,1]))
coord <- as.data.frame(xyFromCell(masked, valid_pixel))

# 3). Compute 10th and 90th percentiles of each day across a 10-day window over 20 years 

# Convert column names into dates
dates <- as.Date(colnames(df_normals))

# Associate each date with the day of a year (ex. 32 -> February 1)
day_year <- yday(dates) ### vector of days
day_map <- tibble(date = dates, day = day_year) ### associated data frame

# Percentiles for each day of the year
percentiles_normals <- function(target_day, window = 5) {
  # Compute 10-day window
  window_day <- ((target_day - window):(target_day + window) - 1) %% 366 + 1
  
  # Select dates corresponding to the window
  selected_dates <- day_map |>
    filter(day_year %in% window_day) |>
    pull(date) |>
    as.character()
  
  selected_data <- df_normals[, selected_dates, drop = FALSE]
  
  # Compute percentiles for the window of the target day
  apply(selected_data, 1, function(x) quantile(x, probs = c(0.1, 0.9), na.rm = TRUE))
}

# Compute for all days of a year
day_normals <- lapply(1:366, percentiles_normals)

# 4). Create data frame containing the percentiles

# Transpose and assemble list day_normals to create a 1142 x (2*366) matrix 
percentiles_mat <- purrr::map(day_normals, ~ t(.x))  

# Extract each day's percentiles
p10_list <- purrr::map(percentiles_mat, ~ .x[,1])
p90_list <- purrr::map(percentiles_mat, ~ .x[,2])

# Convert lists to data frames
p10_df <- do.call(cbind, p10_list)
p90_df <- do.call(cbind, p90_list)

# Rename data frame colums
colnames(p10_df) <- paste0("p10_day_", 1:366)
colnames(p90_df) <- paste0("p90_day_", 1:366)

# Bind the two data frames to obtain the final data frame with the season normals
df_normals_final <- cbind(coord, p10_df, p90_df)

saveRDS(df_normals_final, "data/prepared/normals.rds")

#===============================================================================
# 4) Create warm spell and cold spell duration indexes ------
#===============================================================================

warm_spells_indices()

cold_spells_indices()




