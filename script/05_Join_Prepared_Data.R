#===============================================================================
# Description: Joining together the previously prepared monthly datasets into one 
# large dataframe with communes in rows and monthly aggregates in columns
#===============================================================================

# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ecmwfr, tidyverse, dplyr, terra, maps, here, ncdf4, raster, climate, devtools, lubridate,
               sf, sp, rnaturalearth, Matrix, readr, zoo, gstat, purrr)

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
# 1)  Prepare the data sets for joining ------
#===============================================================================

# Renaming columns of each dataset according to the month covered
rename_columns <- function(source_dir, prefix = TRUE) {
  
  # 1). Find all .rds file for a given month and load them 
  files <- list.files(source_dir, full.names = TRUE, pattern = "\\.rds$")
  files <- files[basename(files) != "normals.rds"]
  
  # 2). Rename each column of each file according to the year-month of the file
  for (file in files) {
    # Get file name without extension
    base_name <- tools::file_path_sans_ext(basename(file))
    
    # Extract year-month from filename (for example "2023-01")
    month_match <- regmatches(base_name, regexpr("\\d{4}-\\d{2}", base_name))
    
    df <- readRDS(file)
    
    # Separate columns to rename and not to rename
    static_cols <- c("x", "y", "ID", "nom", "geometry", "file")
    to_rename <- setdiff(names(df), static_cols)
    
    new_names <- if (prefix) {
      paste(month_match, to_rename, sep = "_")
    } else {
      paste(to_rename, month_match, sep = "_")
    }
    
    names(df)[match(to_rename, names(df))] <- new_names
    
    # Save over original file
    saveRDS(df, file)
  }
}

rename_columns(dir$prepared)

#===============================================================================
# 2)  Join the data sets together ------
#===============================================================================

library(dplyr)
library(purrr)

join_dataset <- function(source_dir, by_cols = c("nom", "x", "y")) {
  files <- list.files(source_dir, full.names = TRUE, pattern = "\\.rds$")
  files <- files[basename(files) != "normals.rds"]
  
  df <- lapply(files, readRDS)
  
  base <- df[[1]]
  
  for (i in 2:length(df)) {
    next_df <- df[[i]]
    
    # Drop geometry from non-base to avoid multiple geometry columns
    if ("geometry" %in% names(next_df)) {
      next_df <- st_drop_geometry(next_df)
    }
  }
  
  final_df <-  left_join(base, next_df, by = by_cols)
  
  return(final_df)
}

final_dataset <- join_dataset(dir$prepared)

saveRDS(final_dataset, "prepared/final_commune_era5_dataset.rds")

