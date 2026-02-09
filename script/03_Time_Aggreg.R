#===============================================================================
# Description: Create yearly, quarterly and monthly weather data
# And combine with climate normals from the chosen reference period
# How to use this script:
# (i) In step 1), specify your variable, your period of data, and your reference period for climate normals
# Note that climat normals files for that periods must have been previously generated using script "02_WeatherNormals.R"
# (ii) Then, all steps of the script can be executed, and data will be saved in dir$prepared
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

#===============================================================================
# 1) First precise what weather statistic you want, specify your
# period of weather analysis, and specifying your climate historical reference period------
#===============================================================================

params <- tribble( # Comment to remove statistics you don't want
  ~variable,        ~statistic,
  "2m_temperature", "daily_mean",
  # "2m_temperature", "daily_minimum",
  # "2m_temperature", "daily_maximum",
  "total_precipitation", "daily_sum"
)

# Period of analysis
start_date <- ymd("2015-01-01")
end_date <- ymd("2023-12-31")

# Reference period:
# Must be sure that there are data files for climate normales on this period
# If not, generate the normales with script "02_WeatherNormals.R"
start_date_normals <- "1971"
end_date_normals <- "2000"

#===============================================================================
# 2) Combine all ncdf files in a dataframe to get all daily observations 
# in a single object and then compute monthly, quarterly and yearly aggregates
# for temperature and precipitation ------
#===============================================================================

# Loop for each row in params
weather_list <- params %>%
  filter(variable %in% c("2m_temperature")) %>% 
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
    
    })
  )

temperature_df <- weather_list %>% 
  filter(variable == "2m_temperature") %>% 
  dplyr::select(output) %>%
  unnest(cols = output) %>% 
  mutate(
    year = year(date),
    quarter = quarter(date),
    month = month(date),
    day = day(date),
    doy = yday(date),
    temp_C = `2m_temperature_daily_mean` - 273.15
  )

precipitation_df <- weather_list %>% 
  filter(variable == "total_precipitation") %>% 
  dplyr::select(output) %>%
  unnest(cols= output) %>% 
  mutate(
    year = year(date),
    quarter = quarter(date),
    month = month(date),
    day = day(date),
    doy = yday(date)
  )
   
#===============================================================================
# 3) Load climatic normales data files  ------
#===============================================================================

# Find all prepared files
prepared_files <- list.files(here(dir$prepared))
# Filter to retain files with "normales" in their name
normals_files <- prepared_files[grepl("climate_normales", prepared_files)]
# Filter to retain files with reference period corresponding to the one specified
normals_ref_files <- normals_files[grepl(paste0(start_date_normals, "_", end_date_normals), normals_files)]

# Load temperature normals
temperature_normals_percentiles <- readRDS(here(dir$prepared, normals_ref_files[grepl("percentiles", normals_ref_files)]))
temperature_normals <- readRDS(here(dir$prepared, normals_ref_files[grepl("temperature", normals_ref_files) & !( grepl("percentile", normals_ref_files))]))

# Load precipitation normals
precipitation_normals <- readRDS(here(dir$prepared, normals_ref_files[grepl("precipitation", normals_ref_files)]))

# Add rolling percentiles to temperature df
temperature_df <- temperature_df %>%
  left_join(temperature_normals_percentiles, by = c("x", "y", "month", "day", "variable", "statistic")) 

# Add percentiles to precipitation df
precipitation_normals_percentiles <- precipitation_normals %>%
  group_by(x,y) %>% 
  slice(1) %>% 
  dplyr::select(x, y, PREC0.90P, PREC0.95P, PREC0.99P)
precipitation_df <- precipitation_df %>%
  left_join(precipitation_normals_percentiles, by = c("x", "y"))

#===============================================================================
# 4) Create yearly, monthly and quarterly temperature data ------
#===============================================================================

df <- temperature_df 
# Transform in data table
df <-  setDT(df)

# Monthly stats (avec data.table)
monthly <- df[
  ,
  .(
    min = min(temp_C, na.rm = TRUE),
    max = max(temp_C, na.rm = TRUE),
    mean = mean(temp_C, na.rm = TRUE),
    sd = sd(temp_C, na.rm = TRUE),
    cold_days = sum(temp_C < p0.10, na.rm=T), # remove NA as there is one for leap years (29th of feb has no percentiles)
    hot_days = sum(temp_C > p0.90, na.rm=T)
  ),
  by = .(x, y, year, month)
][, freq := "monthly"]

# Quartely stats
quarterly <- df[
  ,
  .(
    min = min(temp_C, na.rm = TRUE),
    max = max(temp_C, na.rm = TRUE),
    mean = mean(temp_C, na.rm = TRUE),
    sd = sd(temp_C, na.rm = TRUE),
    cold_days = sum(temp_C < p0.10, na.rm = T), # remove NA as there is one for leap years (29th of feb has no percentiles)
    hot_days = sum(temp_C > p0.90, na.rm = T)
  ),
  by = .(x, y, year, quarter)
][, freq := "quarterly"]

# Yearly stats
yearly <- df[
  ,
  .(
    min = min(temp_C, na.rm = TRUE),
    max = max(temp_C, na.rm = TRUE),
    mean = mean(temp_C, na.rm = TRUE),
    sd = sd(temp_C, na.rm = TRUE),
    cold_days = sum(temp_C < p0.10, na.rm=T), # remove NA as there is one for leap years (29th of feb has no percentiles)
    hot_days = sum(temp_C > p0.90, na.rm=T)
  ),
  by = .(x, y, year)
][, freq := "yearly"]


# Compute number of days within temperature bins
# Define bins: cut() returns factor, so force levels
bins <- c(-100, 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, Inf)
bin_labels <- c("lower_0", "0_3", "3_6","6_9", "9_12", "12_15", "15_18", "18_21",
                "21_24", "24_27", "27_30", "greater_30")

df[, temp_bin := cut(
  temp_C,
  breaks = bins,
  labels = bin_labels,
  right = FALSE, 
  include.lowest = TRUE
)]

# Helper: aggregate, pivot wider
bin_counts <- function(by_cols, freq_name) {
  # by_cols = "month"
  # freq_name = "monthly"
  tmp <- df[
    ,
    .(days = .N),
    by = c("x", "y", "year", by_cols, "temp_bin")
  ]
  
  
  tmp2 <- dcast(tmp, ... ~ temp_bin, value.var = "days", fill = 0)
  tmp2[, freq := freq_name]
  
  setnames(tmp2, old = bin_labels, new = paste0("days_bin_", bin_labels))
  
  tmp2
}

bins_month <- bin_counts("month", "monthly")
bins_quarter <- bin_counts("quarter", "quarterly")
bins_year <- bin_counts("year", "yearly")

# Add bins data to monthly, quarterly and yearly data
monthly <- monthly %>%
  left_join(bins_month, by = c("x", "y", "year", "month", "freq"))
quarterly <- quarterly %>%
  left_join(bins_quarter, by = c("x", "y", "year", "quarter", "freq"))
yearly <- yearly %>% 
  left_join(bins_year, by = c("x", "y", "year", "freq")) %>% 
  dplyr::select(- year.1)

# Add climate normals to the table:
temperature_normals_renamed <- temperature_normals %>%
  rename_with(
    ~ paste0("reference_", .),
    -c(x, y, freq, quarter, month, reference_period, statistic, variable)
  )
# For monthly data
temperature_normals_monthly <- temperature_normals_renamed %>% 
  filter(freq == "monthly") %>% 
  dplyr::select(-quarter)
monthly_weather_df <- monthly %>% 
  left_join(temperature_normals_monthly, by = c("x", "y", "month", "freq"))

# For quarterly data
temperature_normals_quarterly <- temperature_normals_renamed %>% 
  filter(freq == "quarterly") %>% 
  dplyr::select(-month)
quarterly_weather_df <- quarterly %>% 
  left_join(temperature_normals_quarterly, by = c("x", "y", "quarter", "freq"))

# For yearly
temperature_normals_yearly <- temperature_normals_renamed %>% 
  filter(freq == "historical") %>% 
  mutate(freq = "yearly") %>% 
  dplyr::select(-c(quarter, month))
yearly_weather_df <- yearly %>%
  left_join(temperature_normals_yearly, by = c("x", "y", "freq"))

# Combine results
temperature_df <- bind_rows(
  monthly_weather_df,
  quarterly_weather_df,
  yearly_weather_df
  ) %>%
  arrange(freq, year, quarter, month, x,y) %>% 
  dplyr::select(x, y, year, month, quarter, freq, variable, statistic, reference_period, everything())

#===============================================================================
# 4) Create yearly, monthly and quarterly precipitation data ------
#===============================================================================

# Transform in data table
df_precipitation <-  setDT(precipitation_df)

# Create binary variable for dry days, wet days (PRCP > 0.001m (1mm)), very wet and very very wet days
df_precipitation[, dry_day := ifelse(total_precipitation_daily_sum < 0.001, TRUE, FALSE)]
df_precipitation[, wet_day := ifelse(total_precipitation_daily_sum >= 0.001, TRUE, FALSE)]
df_precipitation[, wet10mm := total_precipitation_daily_sum >= 0.010]
df_precipitation[, wet20mm := total_precipitation_daily_sum >= 0.020]

# Compute number of consecutive dry days within temperature bins
# First arrange data by (x,y) and then date
setorder(df_precipitation, x, y, date)

# Run ID over the entire series for each cell
df_precipitation[, dry_run_id := rleid(x, y, dry_day)]
df_precipitation[, wet_run_id := rleid(x, y, wet_day)]

# Dry spell length within each cell by month, quarter and year
df_precipitation[, dry_spell_in_month := seq_len(.N) * dry_day, by = .(x, y, year, month, dry_run_id)]
df_precipitation[, dry_spell_in_quarter := seq_len(.N) * dry_day, by = .(x, y, year, quarter, dry_run_id)]
df_precipitation[, dry_spell_in_year := seq_len(.N) * dry_day, by = .(x, y, year, dry_run_id)]

# Wet spell length within each cell by month, quarter and year
df_precipitation[, wet_spell_in_month := seq_len(.N) * wet_day, by = .(x, y, year, month, dry_run_id)]
df_precipitation[, wet_spell_in_quarter := seq_len(.N) * wet_day, by = .(x, y, year, quarter, dry_run_id)]
df_precipitation[, wet_spell_in_year := seq_len(.N) * wet_day, by = .(x, y, year, dry_run_id)]


# Monthly stats
monthly <- df_precipitation[
  ,
  .(
    wet_days = sum(wet_day, na.rm = TRUE),
    dry_days = sum(dry_day, na.rm = TRUE),
    total_precip_wet_days = sum(fifelse(wet_day, total_precipitation_daily_sum, 0), na.rm = TRUE),
    wet10mm_days = sum(wet10mm, na.rm = TRUE),
    wet20mm_days = sum(wet20mm, na.rm = TRUE),
    R90p =  sum(total_precipitation_daily_sum[total_precipitation_daily_sum > PREC0.90P], na.rm = TRUE),
    R95p =  sum(total_precipitation_daily_sum[total_precipitation_daily_sum > PREC0.95P], na.rm = TRUE),
    R99p =  sum(total_precipitation_daily_sum[total_precipitation_daily_sum > PREC0.99P], na.rm = TRUE),
    Max_CDD = max(dry_spell_in_month, na.rm = T),
    Max_CWD = max(wet_spell_in_month, na.rm = T)
  ),
  by = .(x, y, year, month)
][, freq := "monthly"]

# Monthly stats
quarterly <- df_precipitation[
  ,
  .(
    wet_days = sum(wet_day, na.rm = TRUE),
    dry_days = sum(dry_day, na.rm = TRUE),
    total_precip_wet_days = sum(fifelse(wet_day, total_precipitation_daily_sum, 0), na.rm = TRUE),
    wet10mm_days = sum(wet10mm, na.rm = TRUE),
    wet20mm_days = sum(wet20mm, na.rm = TRUE),
    R90p =  sum(total_precipitation_daily_sum[total_precipitation_daily_sum > PREC0.90P], na.rm = TRUE),
    R95p =  sum(total_precipitation_daily_sum[total_precipitation_daily_sum > PREC0.95P], na.rm = TRUE),
    R99p =  sum(total_precipitation_daily_sum[total_precipitation_daily_sum > PREC0.99P], na.rm = TRUE),
    Max_CDD = max(dry_spell_in_quarter, na.rm = T),
    Max_CWD = max(wet_spell_in_quarter, na.rm = T)
  ),
  by = .(x, y, year, quarter)
][, freq := "quarterly"]

# Yearly stats
yearly <- df_precipitation[
  ,
  .(
    wet_days = sum(wet_day, na.rm = TRUE),
    dry_days = sum(dry_day, na.rm = TRUE),
    total_precip_wet_days = sum(fifelse(wet_day, total_precipitation_daily_sum, 0), na.rm = TRUE),
    wet10mm_days = sum(wet10mm, na.rm = TRUE),
    wet20mm_days = sum(wet20mm, na.rm = TRUE),
    R90p =  sum(total_precipitation_daily_sum[total_precipitation_daily_sum > PREC0.90P], na.rm = TRUE),
    R95p =  sum(total_precipitation_daily_sum[total_precipitation_daily_sum > PREC0.95P], na.rm = TRUE),
    R99p =  sum(total_precipitation_daily_sum[total_precipitation_daily_sum > PREC0.99P], na.rm = TRUE),
    Max_CDD = max(dry_spell_in_year, na.rm = T),
    Max_CWD = max(wet_spell_in_year, na.rm = T)
  ),
  by = .(x, y, year)
][, freq := "yearly"]


# Add precipitation normales to the table:
precipitation_normals_renamed <- precipitation_normals %>%
  rename_with(
    ~ paste0("reference_", .),
    -c(x, y, freq, quarter, month, reference_period, period)
  )
# For monthly data
precipitation_normals_monthly <- precipitation_normals_renamed %>% 
  filter(freq == "monthly") %>% 
  dplyr::select(-quarter)
monthly_precipitation_df <- monthly %>% 
  left_join(precipitation_normals_monthly, by = c("x", "y", "month", "freq"))

# For quarterly data
precipitation_normals_quarterly <- precipitation_normals_renamed %>% 
  filter(freq == "quarterly") %>% 
  dplyr::select(-month)
quarterly_precipitation_df <- quarterly %>% 
  left_join(precipitation_normals_quarterly, by = c("x", "y", "quarter", "freq"))

# For yearly
precipitation_normals_yearly <- precipitation_normals_renamed %>% 
  filter(freq == "yearly") %>%   
  dplyr::select(-c(quarter, month))
yearly_precipitation_df <- yearly %>%
  left_join(precipitation_normals_yearly, by = c("x", "y", "freq"))

# Combine results
precipitation_df <- bind_rows(
  monthly_precipitation_df,
  quarterly_precipitation_df,
  yearly_precipitation_df
) %>%
  mutate(
    variable = "total_precipitation",
    statistic = "daily_sum"
  ) %>% 
  arrange(freq, year, quarter, month, x,y) %>% 
  dplyr::select(x, y, year, month, quarter, freq, variable, statistic, reference_period, everything())


#===============================================================================
# 5) Save data ------
#===============================================================================

saveRDS(temperature_df, 
       file = here(
         dir$prepared, 
         paste("weather_temperature_",
               year(start_date),
               "_",
               year(end_date),
               "_reference_",
               start_date_normals,
               "_",
               end_date_normals,
               ".rds",
               sep="")
         )
       )

saveRDS(precipitation_df, 
        file = here(
          dir$prepared, 
          paste("weather_precipitation_",
                year(start_date),
                "_",
                year(end_date),
                "_reference_",
                start_date_normals,
                "_",
                end_date_normals,
                ".rds",
                sep="")
        )
)
