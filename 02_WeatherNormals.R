#===============================================================================
# Description: Create daily weather normals for each grid cells
#===============================================================================

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

params <- tribble(
  ~variable,        ~statistic,
  "2m_temperature", "daily_mean",
  "2m_temperature", "daily_minimum",
  "2m_temperature", "daily_maximum",
  "total_precipitation", "daily_sum"
)

#===============================================================================
# 1) Combine all ncdf files in a dataframe to get all daily observations 
# in a single object and then compute historical weather for temperature ------
#===============================================================================

# Temperature first
start_date <- ymd("2020-01-01")
end_date <- ymd("2023-12-31")

# Loop for each row in params
climate_normals_temperature <- params %>%
  filter(
    variable == "2m_temperature"
  ) %>% 
  mutate(
    output = map2(variable, statistic, ~{
      
      # 1) Process ERA5 files
      df_weather <- process_era5_files(
        source_dir = dir$source,
        variable = .x,
        statistic = .y,
        start_date = start_date,
        end_date = end_date
      ) %>%
        do.call(rbind, .) %>%
        mutate(variable = .x, statistic = .y)
      
      # 2) Compute normals function in script 00_Functions.R
      compute_temperature_normals_dt(
        df = df_weather,
        value_col = paste(.x, .y, sep = "_"),
        rolling_window = 7,
        percentiles = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
      )
    })
  )


# Combine results
# combine all rolling tables, keeping .id = index
all_rolling <- map_dfr(climate_normals_temperature$output, "rolling", .id = "id")

# join with params by id
# First: add id to params too
params2 <- params %>%
  filter(
    variable == "2m_temperature"
  ) %>% 
  mutate(id = as.character(row_number()))

# Now left join!
rolling_temperature_percentiles_df <- all_rolling %>%
  left_join(params2, by = "id") %>% 
  rename(month = target_month,
         day = target_day) %>%
  rename_with(
    .fn = function(x) {
      percent_cols <- gsub("%", "", x)
      paste0("p0.", sprintf("%02d", as.numeric(percent_cols)))
    },
    .cols = matches("^[0-9]+%$")
  ) %>% 
  select(-id) %>% 
  arrange(day, month, x, y, statistic)

monthly_normales_temperature_df <- map_dfr(climate_normals_temperature$output, "monthly", .id = "id") %>%
  left_join(params2, by = "id") %>% 
  as.tibble() %>% 
  mutate(
    min_monthly = min,
    max_monthly = max,
    mean_monthly = mean,
    sd_monthly = sd
  ) %>% 
  dplyr::select(-id, -min, -max, -mean, -sd) %>% 
  arrange(month, x, y, statistic)

#===============================================================================
# 2) Combine all ncdf files in a dataframe to get all daily observations 
# in a single object and then compute historical weather for precipitations ------
#===============================================================================

# Precipitations second
start_date <- ymd("2020-01-01")
end_date <- ymd("2023-12-31")

# Loop for each row in params
climate_normals_temperature <- params %>%
  filter(
    variable == "total_precipitation"
  ) %>% 
  mutate(
    output = map2(variable, statistic, ~{
      
      # 1) Process ERA5 files
      df_weather <- process_era5_files(
        source_dir = dir$source,
        variable = .x,
        statistic = .y,
        start_date = start_date,
        end_date = end_date
      ) %>%
        do.call(rbind, .) %>%
        mutate(variable = .x, statistic = .y)
      
      # 2) Compute normals using function in script 00_Functions.R
      compute_precipitations_normals_dt(
        df = df_weather,
        value_col = paste(.x, .y, sep = "_"),
        rolling_window = 7,
        percentiles = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
      )
    })
  )

#===============================================================================
# 3) Save historical weather ------
#===============================================================================

saveRDS(
  rolling_temperature_percentiles_df,
  file = here(dir$prepared, "climate_normals_temperature_percentiles.rds")
)
saveRDS(
  monthly_normales_temperature_df,
  file = here(dir$prepared, "monthly_normales_temperature.rds")
)
