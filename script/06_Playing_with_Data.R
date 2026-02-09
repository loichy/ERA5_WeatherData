#======================================================================================
# Using scripts from 03_Time_Aggreg.R to play with the temperature data, make maps, ...
#======================================================================================

# Clean memory 
rm(list=ls())
gc()

# Load package
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ecmwfr, tidyverse, dplyr, purrr, terra, maps, here, ncdf4, raster, climate, devtools, lubridate,
               sf, sp, rnaturalearth, rnaturalearthdata, Matrix, slider, data.table)

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

#==========================================================
# We use the dataframe we saved in script 03_Time_Aggreg.R
#==========================================================

df <- readRDS("data/prepared/weather_temperature_2015_2023_reference_1971_2000.rds")


# Let's see how temperatures evolve in Mendive (Basque Country) 
# Mendive : x = -1 & y = 43

mendive_temp <- df %>%
  filter(x == -1, y == 43,
         freq == "quarterly") %>%
  group_by(year, quarter)

ggplot(mendive_temp, aes(x = year, y = mean, color = as.character(quarter))) +
  geom_point() + geom_line() + theme_minimal() +
  labs(x = "Year",
       y = "Mean temperature", 
       title = "Evolution of the mean temperature in Mendive by quarter",
       color = "Quarter")
