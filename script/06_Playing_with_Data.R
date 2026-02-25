#======================================================================================
# Using script from 03_Time_Aggreg.R to play with the temperature data, make maps, ...
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
normales_temperature_df <- readRDS("data/prepared/climate_normales_temperature_1971_2000.rds")
daily_temp_df <- readRDS("data/prepared/daily_temperature_2015_2023_reference_1971_2000.rds")

select <- dplyr::select 


# Let's see how temperatures evolve in Mendive (Basque Country) 
# Mendive : x = -1 & y = 43

mendive_temp <- df %>%
  filter(x == -1, y == 43,
         freq == "quarterly", 
         ! is.na(quarter)) %>%
  select(year, min, max, mean, reference_mean, reference_max, reference_min, quarter)

ggplot(mendive_temp, aes(x = year)) +
  geom_line(aes(y = min, color = "Minimum")) + 
  geom_line(aes(y = reference_min, color = "Minimum reference")) +
  scale_color_manual(values = c("blue", "skyblue")) +
  facet_wrap(~ quarter) +
  theme_minimal() +
  labs(title = "Evolution of minimal temperatures by quarter in Mendive",
       x = "Year",
       y = "Temperature in C°",
       color = "Color")

ggplot(mendive_temp, aes(x = year)) +
  geom_line(aes(y = max, color = "Maximum")) +
  geom_line(aes(y = reference_max, color = "Maximum reference")) +
  scale_color_manual(values = c("darkorange1", "red2")) +
  facet_wrap(~ quarter) +
  theme_minimal() +
  labs(title = "Evolution of maximal temperatures by quarter in Mendive",
       x = "Year",
       y = "Temperature in C°", 
       color = "Color")

ggplot(mendive_temp, aes(x = year)) +
  geom_line(aes(y = mean, color = "Mean")) +
  geom_line(aes(y = reference_mean, color = "Mean reference")) +
  scale_color_manual(values = c("chartreuse", "forestgreen")) +
  facet_wrap(~ quarter) +
  theme_minimal() +
  labs(title = "Evolution of mean temperatures by quarter in Mendive",
       x = "Year",
       y = "Temperature in C°",
       color = "Color")

# Let's see how the number of hot days has evolved from 2015 to 2023 for France


ggplot(df %>% filter(year %in% c(2015, 2023)),
  aes(x = x, y = y, color = hot_days, fill = hot_days)) +
  geom_tile() + 
  scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(~ year) +  
  theme_void() +
  labs(title = "Comparison of Hot Days",
    subtitle = "Years: 2015 vs 2023",
    fill = "Hot days",
    color = "Hot days")

ggplot(df, aes(x = x, y = y, color = hot_days, fill = hot_days)) +
  geom_tile() + 
  facet_wrap(~ year) +
  scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  labs(title = "Number of hot days by year",
       fill = "Hot days",
       color = "Hot days")

# Which seasons are the most impacted?

ggplot(df %>%
         select(x, y, year, quarter, hot_days) %>%
         filter(year == 2015, ! is.na(quarter)), 
       aes(x = x, y = y, color = hot_days, fill = hot_days)) +
  geom_tile() + 
  facet_wrap( ~ quarter) + 
  scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  labs(title = "Number of hot days in 2015 by quarter",
       fill = "Hot days",
       color = "Hot days")

ggplot(df %>%
         select(x, y, year, quarter, hot_days) %>%
         filter(year == 2023, ! is.na(quarter)), 
       aes(x = x, y = y, color = hot_days, fill = hot_days)) +
  geom_tile() + 
  facet_wrap( ~ quarter) + 
  scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  labs(title = "Number of hot days in 2023 by quarter",
       fill = "Hot days",
       color = "Hot days")

# Which months are the most impacted?

ggplot(df %>%
         select(x, y, year, month, hot_days) %>%
         filter(year == 2015, ! is.na(month)), 
       aes(x = x, y = y, color = hot_days, fill = hot_days)) +
  geom_tile() + 
  facet_wrap( ~ month) + 
  scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  labs(title = "Number of hot days in 2015 by month",
       fill = "Hot days",
       color = "Hot days")

ggplot(df %>%
         select(x, y, year, month, hot_days) %>%
         filter(year == 2023, ! is.na(month)), 
       aes(x = x, y = y, color = hot_days, fill = hot_days)) +
  geom_tile() + 
  facet_wrap( ~ month) + 
  scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  labs(title = "Number of hot days in 2023 by month",
       fill = "Hot days",
       color = "Hot days")


# ==============================================
# Let's create a few more indices (ETCCDI)
# ==============================================

# 1- We compute GSL : the growing season length (n°5)

GSL_daily <- daily_temp_df %>%
  filter(x == 6.75, y == 45.25)


calc_gsl_fast <- function(temp, doy) {
  # Let's make sure the data is in the right order
  # (Even though is shoud already be thanks to the doy)
  ord <- order(doy)
  t <- temp[ord]
  d <- doy[ord]
  
  # --- FIND THE START ---
  # We are looking for 6 consecutive days > 5°C
  runs_start <- rle(t > 5)
  start_indices <- cumsum(runs_start$lengths) - runs_start$lengths + 1
  
  start_doy <- NA
  for(i in seq_along(runs_start$values)) {
    if(runs_start$values[i] == TRUE && runs_start$lengths[i] >= 6) {
      start_doy <- d[start_indices[i]]
      break
    }
  }
  
  # If there is no growing season starting, we return 0
  if(is.na(start_doy)) return(0)
  
  # --- FIND THE END ---
  # We can oly look at temperatures after the 1st of July (doy 182)
  post_july_mask <- d >= 182
  t_post <- t[post_july_mask]
  d_post <- d[post_july_mask]
  
  end_doy <- max(d) # Default value if winter doesn't come back before 31/12
  
  runs_end <- rle(t_post < 5)
  end_indices <- cumsum(runs_end$lengths) - runs_end$lengths + 1
  
  for(i in seq_along(runs_end$values)) {
    if(runs_end$values[i] == TRUE && runs_end$lengths[i] >= 6) {
      end_doy <- d_post[end_indices[i]]
      break
    }
  }
  
  return(end_doy - start_doy)
}

# GSL by location and year
gsl_results <- daily_temp_df %>%
  group_by(x, y, year) %>%
  summarize(gsl_days = calc_gsl_fast(temp_C, doy), .groups = "drop")


# Il existe aussi un package : ClimInd et une fonction gsl 
# pour faire ça mais j'arrive pas encore à le faire marcher
# library(ClimInd)

# Time to visualize the results !


ggplot(gsl_results %>% filter(year %in% c(2015, 2023)),
       aes(x = x, y = y, color = gsl_days, fill = gsl_days)) +
  geom_tile() + 
  scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(~ year) +  
  theme_void() +
  labs(title = "Comparison of Growing season length",
       subtitle = "Years: 2015 vs 2023",
       fill = "Growing season days",
       color = "Growing season days")

gsl_diff <- gsl_results %>% # To compute the difference between 2023/2015
  filter(year %in% c(2015, 2023)) %>%
  pivot_wider(names_from = year,
              values_from = gsl_days,
              names_prefix = "gsl_") %>%
  mutate(dif = gsl_2023 - gsl_2015)

ggplot(gsl_diff,
       aes(x = x, y = y, color = dif, fill = dif)) +
  geom_tile() + 
  scale_fill_viridis_c(option = "turbo") +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  labs(title = "Evolution of Growing season length",
       subtitle = "Years: 2015 to 2023",
       fill = "Growing season difference",
       color = "Growing season difference")

summary(gsl_diff$dif) # Le représenter dans un joli tableau

