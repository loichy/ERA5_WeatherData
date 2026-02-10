#===============================================================================
# Description: Create daily weather normals for each grid cells
#===============================================================================

# Clean memory 
rm(list=ls())
gc()

# Load package
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ecmwfr, tidyverse, dplyr, purrr, terra, maps, here, ncdf4, raster, climate, devtools, lubridate,
               sf, sp, rnaturalearth, rnaturalearthdata, Matrix, slider, data.table)#test
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
start_date <- ymd("1971-01-01")
end_date <- ymd("2000-12-31")

# Loop for each row in params
climate_normals_temperature <- params %>%
  filter(
    variable == "2m_temperature",
    statistic == "daily_mean"
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
all_normals <- map_dfr(climate_normals_temperature$output, "normals", .id = "id") 

# join with params by id
# First: add id to params too
params2 <- params %>%
  filter(
    variable == "2m_temperature",
    statistic == "daily_mean"
  ) %>% 
  mutate(id = as.character(row_number()))

# Now left join:
normales_temperature_df <- all_normals %>% 
  left_join(params2, by = "id") %>% 
  dplyr::select(-id) %>% 
  arrange(freq, quarter, month, x, y) %>% 
  mutate(reference_period = paste(year(start_date), year(end_date), sep = "_"))

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
  dplyr::select(-id) %>% 
  arrange(day, month, x, y, statistic) %>% 
  mutate(reference_period = paste(year(start_date), year(end_date), sep = "_"))


# check_temp <- nomales_temperature_df %>%
#   filter(freq == "historical") %>% 
#   mutate(
#    hot_days = days_bin_24_27 + days_bin_27_30 + days_bin_greater_30 
#   )
#  
# ggplot(data = check_temp, aes(x=x, y=y, color = hot_days, fill = hot_days))+
#   geom_tile() +
#   scale_fill_viridis_c() +
#   scale_color_viridis_c() +
#   theme_bw() +
#   labs(x = "Longitude",
#        y = "Latitude") +
#   theme(legend.position = "bottom")


#===============================================================================
# 2) Combine all ncdf files in a dataframe to get all daily observations 
# in a single object and then compute historical weather for precipitations ------
#===============================================================================

# Precipitations second
start_date <- ymd("1971-01-01")
end_date <- ymd("2000-12-31")

# Loop for each row in params
climate_normals_precipitation <- params %>%
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
        mutate(variable = .x, 
               statistic = .y)

      # 2) Compute normals using function in script 00_Functions.R
      compute_precipitation_normals_dt(
        df = df_weather,
        value_col = paste(.x, .y, sep = "_")
      )
    })
  )


normales_precipitation_df <- map_dfr(climate_normals_precipitation$output, "result", .id = "id") %>% 
  mutate(reference_period = paste(year(start_date), year(end_date), sep = "_"))

# check_prec <- normales_precipitation_df %>% 
#   filter(freq == "yearly") 
# 
# ggplot(data = check_prec, aes(x=x, y=y, color = avg_wet_days, fill = avg_wet_days))+
#   geom_tile() +
#   scale_fill_viridis_c() +
#   scale_color_viridis_c() +
#   theme_bw() +
#   labs(x = "Longitude",
#        y = "Latitude") +
#   theme(legend.position = "bottom")

#===============================================================================
# 3) Save historical weather ------
#===============================================================================

saveRDS(
  rolling_temperature_percentiles_df,
  file = here(
    dir$prepared, 
    paste("climate_normales_temperature_percentiles_", 
          year(start_date),
          "_",
          year(end_date),
          ".rds", sep = ""
          )
  )
)
saveRDS(
  normales_temperature_df,
  file = here(
    dir$prepared, 
    paste("climate_normales_temperature_",
          year(start_date),
          "_",
          year(end_date),
          ".rds", sep = ""
          )
    )
)
saveRDS(
  normales_precipitation_df,
  file = here(
    dir$prepared, 
    paste("climate_normales_precipitation_",
          year(start_date),
          "_",
          year(end_date),
          ".rds", sep = ""
          )
    )
)

