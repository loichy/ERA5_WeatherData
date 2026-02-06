#===============================================================================
# Description: Functions to download and aggregate ERA5-LAND data
#===============================================================================

# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ecmwfr, tidyverse, dplyr, terra, maps, here, ncdf4, raster, climate, devtools, lubridate,
               sf, sp, rnaturalearth, Matrix, readr, zoo, gstat, FNN)

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

#===============================================================================

# Download function over multiple months for source data

download_era5_ecmwfr <- function(variable,
                                 start_date,
                                 end_date,
                                 statistic,
                                 area = c(52, -6, 41, 10),
                                 # frequency = "6_hourly",
                                 time_zone = "utc+02:00",
                                 dataset = "derived-era5-single-levels-daily-statistics",
                                 product_type = "reanalysis",
                                 user = "ecmwfr",
                                 key,
                                 output_file = NULL,
                                 output_dir) {
  
  # 1). Set API key
  wf_set_key(key = key, user = user)
  
  # 2). Create subdirectory for the variable and statistic
  variable_dir <- file.path(output_dir, variable, statistic)
  if (!dir.exists(variable_dir)) {
    dir.create(variable_dir, recursive = TRUE)
  }
  
  # 3). Generate date groups by month
  dates <- seq(ymd(start_date), ymd(end_date), by = "day")
  date_df <- tibble(date = dates) |>
    mutate(year = format(date, "%Y"),
           month = format(date, "%m"),
           day = format(date, "%d")) |>
    group_by(year, month)
  
  downloaded_files <- c()
  
  # 4). Loop through each group
  for (grp in group_split(date_df)) {
    year <- unique(grp$year)
    month <- unique(grp$month)
    days <- grp$day
    
    # Monthly file name
    if (is.null(output_file)) {
      output_file_month <- paste0("era5_", variable, "_", year, "-", month, "_", statistic, ".nc")
    } else {
      output_file_month <- gsub("\\.nc$", paste0("_", year, "-", month, ".nc"), output_file)
    }
    
    # Construct API request
    request <- list(
      dataset_short_name = dataset,
      product_type = product_type,
      variable = variable,
      year = year,
      month = month,
      day = days,
      daily_statistic = statistic,
      time_zone = time_zone,
      # frequency = frequency,
      area = area,
      format = "netcdf",
      target = output_file_month
    )
    
    # Launch request
    file_path <- wf_request(
      request = request,
      transfer = TRUE,
      user = user, 
      path = variable_dir
    )
    
    downloaded_files <- c(downloaded_files, file_path)
  }
  
  return(invisible(downloaded_files))
}

#===============================================================================

# Helper function: extract dates from NetCDF file
extract_dates <- function(nc_file) {
  # nc_file <- temp_files[1]
  nc <- nc_open(nc_file)
  # Extract the time variable
  time_vals <- ncvar_get(nc, "valid_time")
  time_units_attr <- ncatt_get(nc, "valid_time", "units")
  
  # Extract origin safely
  if (is.list(time_units_attr)) {
    time_units <- time_units_attr$value
  } else {
    time_units <- time_units_attr
  }
  
  origin <- sub(".*since ", "", time_units)
  dates <- as.Date(time_vals, origin = origin)
  nc_close(nc)
  return(dates)
}  #Add france layers and cropped (removed it unfortunatey)

#===============================================================================

# Automation of downloaded data processing and border delimitation

process_era5_files <- function(source_dir, 
                               variable, 
                               statistic,
                               start_date,
                               end_date) {
  # source_dir <- dir$source
  # variable <- "2m_temperature"
  # statistic <- "daily_mean"
  # start_date <- ymd("2020-01-01") # Reference period for weather normals
  # end_date <- ymd("2020-12-31")
  
  # Convert start_date and end_date to Date objects
  start_date <- ymd(start_date) # Reference period for weather normals
  end_date <- ymd(end_date)
  
  # 1). Find all .nc file for a given month 
  nc_files <- list.files(source_dir, pattern = "\\.nc$", full.names = TRUE, recursive = T)
  
  # Filter files by variable and statistic
  nc_files <- nc_files[grepl(paste0(variable, ".*", statistic), basename(nc_files))]
  
  # Filter files on the reference period
  nc_files <- nc_files %>%
    tibble(path = .) %>%
    mutate(
      # Extract YYYY-MM from filename
      date_str = str_extract(path, "\\d{4}-\\d{2}"),
      date = ymd(paste0(date_str, "-01"))
    )  %>%
    filter(date >= start_date & date <= end_date) %>%
    pull(path)  # Only return the filtered file paths
  
  # 2). Open France borders shapefile
  # Get France borders with sf package
  france <- rnaturalearth::ne_countries(scale = "medium", country = "France", returnclass = "sf")
  
  # Harmonize projection with ncdf files
  rast_ncdf <- terra::rast(nc_files[1])  # Load the first netCDF file to get the CRS
  france_proj <- st_transform(france, crs(rast_ncdf))
  
  
  # 3). Transform into data frame format for each file
  processed_list <- lapply(nc_files, function(nc_file) {
    # nc_file <- nc_files[1]
    
    # Load netCDF file
    rast_data <- terra::rast(nc_file)  # Load raster
    
    # Crop data using France borders
    cropped <- crop(rast_data, vect(france_proj))
    masked <- mask(cropped, vect(france_proj))
    
    # Get dates in the layers of the raster
    dates <- extract_dates(nc_file)
      
    # Optional: if number of layers in raster doesn't match dates, warn
    if (nlyr(rast_data) != length(dates)) {
      warning(paste("Mismatch in file:", nc_file))
    }
      
    # Truncate or repeat dates to match layers, just in case
    names(rast_data) <- as.character(dates[seq_len(min(length(dates), nlyr(rast_data)))])
    layer_names <- names(rast_data)

    # Crop data using France borders
    cropped <- crop(rast_data, vect(france_proj))
    masked <- mask(cropped, vect(france_proj))
    
    # Transform into data frame format
    df <- as.data.frame(masked, xy = TRUE, na.rm = TRUE)
    
    # Add treated filename
    df$file <- basename(nc_file)
    
    # Pivot data longer
    df <- df %>%
      pivot_longer(cols = -c(x, y, file), 
                   names_to = "date", 
                   values_to = paste(variable, statistic, sep = "_")) %>%
      mutate(date = as.Date(date))  # Convert valid_time to Date type
    
    return(df)
  })
  
  # 4). Return all dataframes
  return(processed_list)
  
}


#===============================================================================

# Creation of climate normals by cells
# Temperature:
compute_temperature_normals_dt <- function(df, value_col, rolling_window = 7, percentiles = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
) {
  
  # variable <- "2m_temperature"
  # statistic <- "daily_mean"
  # value_col <- paste(variable, statistic, sep = "_")
  # rolling_window <- 7
  # df <- climate_normals_temperature$output %>%
  #   do.call(rbind, .)
  # percentiles = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
  
  # Convert to data.table in-place
  setDT(df)
  
  df[, `:=`(
    year = year(date),
    quarter = quarter(date),
    month = month(date),
    day = day(date),
    doy = yday(date),
    temp_C = get(value_col) - 273.15
  )]
  
  # Lookup table for rolling window
  expand_doy_window <- function(target_doy, window = rolling_window) {
    ((target_doy + (-window:window) - 1) %% 365) + 1
  }
  
  days_df <- df[, .(month, day, doy)]
  days_df <- unique(days_df, by = "doy")
  
  lookup <- data.table(target_doy = 1:365)[, .(window_doy = expand_doy_window(target_doy)), by = target_doy]
  lookup <- merge(lookup, days_df[, .(target_doy = doy, target_month = month, target_day = day)], by = "target_doy")
  lookup <- merge(lookup, days_df[, .(window_doy = doy, window_month = month, window_day = day)], by = "window_doy")
  
  # Join in data.table
  df_long <- merge(
    df, lookup,
    by.x = c("month", "day"),
    by.y = c("window_month", "window_day"),
    allow.cartesian = TRUE
  )
  
  # Grouped percentiles: blazing fast!
  rolling <- df_long[
    ,
    as.list(quantile(temp_C, probs = percentiles, na.rm = TRUE, type=1)),
    by = .(x, y, target_month, target_day)
  ]
  
  # Monthly stats
  monthly <- df[
    ,
    .(
      min = min(temp_C, na.rm = TRUE),
      max = max(temp_C, na.rm = TRUE),
      mean = mean(temp_C, na.rm = TRUE),
      sd = sd(temp_C, na.rm = TRUE)
    ),
    by = .(x, y, month)
  ][, freq := "monthly"]
  
  # Quarterly stats
  quarterly <- df[
    ,
    .(
      min = min(temp_C, na.rm = TRUE),
      max = max(temp_C, na.rm = TRUE),
      mean = mean(temp_C, na.rm = TRUE),
      sd = sd(temp_C, na.rm = TRUE)
    ),
    by = .(x, y, quarter)
  ][, freq := "quarterly"]
  
  # Historical stats
  historical <- df[
    ,
    .(
      min = min(temp_C, na.rm = TRUE),
      max = max(temp_C, na.rm = TRUE),
      mean = mean(temp_C, na.rm = TRUE),
      sd = sd(temp_C, na.rm = TRUE)
    ),
    by = .(x, y)
  ][, `:=`(freq = "historical")]
  
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
    
    # Now: average across years to get mean days per period-unit
    tmp2 <- tmp[
      ,
      .(mean_days = mean(days, na.rm = TRUE)),
      by = c("x", "y", by_cols, "temp_bin")
    ]
    
    tmp3 <- dcast(tmp2, ... ~ temp_bin, value.var = "mean_days", fill = 0)
    tmp3[, freq := freq_name]
    
    setnames(tmp3, old = bin_labels, new = paste0("days_bin_", bin_labels))
    
    tmp3
  }
  
  bins_month <- bin_counts("month", "monthly")
  bins_quarter <- bin_counts("quarter", "quarterly")
  bins_year <- bin_counts("year", "yearly")
  bins_historical <- bins_year[, .(
    days_bin_lower_0 = mean(days_bin_lower_0),
    days_bin_0_3 = mean(days_bin_0_3),
    days_bin_3_6 = mean(days_bin_3_6),
    days_bin_6_9 = mean(days_bin_6_9),
    days_bin_9_12 = mean(days_bin_9_12),
    days_bin_12_15 = mean(days_bin_12_15),
    days_bin_15_18 = mean(days_bin_15_18),
    days_bin_18_21 = mean(days_bin_18_21),
    days_bin_21_24 = mean(days_bin_21_24),
    days_bin_24_27 = mean(days_bin_24_27),
    days_bin_27_30 = mean(days_bin_27_30),
    days_bin_greater_30 = mean(days_bin_greater_30)
  ), by = .(x, y)][, `:=`(freq = "historical")]
  
  # -------------- Combine stats + bins --------------
  combined <- rbindlist(list(monthly, quarterly, historical), use.names = TRUE, fill = TRUE)
  combined_bins <- rbindlist(list(bins_month, bins_quarter, bins_historical), use.names = TRUE, fill = TRUE)
  
  final <- merge(combined, combined_bins, by = c("x", "y", "freq", "quarter", "month"), all = TRUE)
  
  return(list(
    rolling = rolling,
    normals = final[]
  ))
}

# Precipitations
compute_precipitation_normals_dt <- function(df, value_col
) {
  
  # variable <- "total_precipitation"
  # statistic <- "daily_sum"
  # value_col <- paste(variable, statistic, sep = "_")
  # df <- climate_normals_temperature$output %>%
  #   do.call(rbind, .)
  
  # Convert to data.table in-place
  setDT(df)
  
  df[, `:=`(
    year = year(date),
    quarter = quarter(date),
    month = month(date),
    day = day(date),
    doy = yday(date)
  )]
  
  # Create binary variable for dry days, wet days (PRCP > 0.001m (1mm)), very wet and very very wet days
  df[, dry_day := ifelse(get(value_col) < 0.001, TRUE, FALSE)]
  df[, wet_day := ifelse(get(value_col) >= 0.001, TRUE, FALSE)]
  df[, wet10mm := get(value_col) >= 0.010]
  df[, wet20mm := get(value_col) >= 0.020]
  
  # Compute precipitations quantile
  quantiles <- df[
    wet_day == TRUE, 
    .(
      # 95th and 99th percentiles only on wet days
      PREC0.99P = quantile(get(value_col), prob = 0.99, na.rm = TRUE),
      PREC0.95P = quantile(get(value_col), prob = 0.95, na.rm = TRUE),
      PREC0.90P = quantile(get(value_col), prob = 0.90, na.rm = TRUE)
  ),
  by = .(x, y)
  ]
  
  # Monthly stats
  by_year_month <- df[
    ,
    .(
      wet_days = sum(wet_day, na.rm = TRUE),
      dry_days = sum(dry_day, na.rm = TRUE),
      total_precip_wet_days = sum(fifelse(wet_day, get(value_col), 0), na.rm = TRUE),
      wet10mm_days = sum(wet10mm, na.rm = TRUE),
      wet20mm_days = sum(wet20mm, na.rm = TRUE)
    ),
    by = .(x, y, year, month)
  ]
  
  monthly <- by_year_month[
    ,
    .(
      avg_wet_days = mean(wet_days, na.rm = TRUE),
      avg_dry_days = mean(dry_days, na.rm = TRUE),
      avg_total_precip_wet_days = mean(total_precip_wet_days, na.rm = TRUE),
      avg_wet10mm_days = mean(wet10mm_days, na.rm = TRUE),
      avg_wet20mm_days = mean(wet20mm_days, na.rm = TRUE)
    ),
    by = .(x, y, month)
  ][, freq := "monthly"]
  
  # Quarterly stats
  by_year_quarter <- df[
    ,
    .(
      wet_days = sum(wet_day, na.rm = TRUE),
      dry_days = sum(dry_day, na.rm = TRUE),
      total_precip_wet_days = sum(fifelse(wet_day, get(value_col), 0), na.rm = TRUE),
      wet10mm_days = sum(wet10mm, na.rm = TRUE),
      wet20mm_days = sum(wet20mm, na.rm = TRUE)
    ),
    by = .(x, y, year, quarter)
  ]
  
  quarterly <- by_year_quarter[
    ,
    .(
      avg_wet_days = mean(wet_days, na.rm = TRUE),
      avg_dry_days = mean(dry_days, na.rm = TRUE),
      avg_total_precip_wet_days = mean(total_precip_wet_days, na.rm = TRUE),
      avg_wet10mm_days = mean(wet10mm_days, na.rm = TRUE),
      avg_wet20mm_days = mean(wet20mm_days, na.rm = TRUE)
    ),
    by = .(x, y, quarter)
  ][, freq := "quarterly"]
  
  # Yearly stat
  by_year <- df[
    ,
    .(
      wet_days = sum(wet_day, na.rm = TRUE),
      dry_days = sum(dry_day, na.rm = TRUE),
      total_precip_wet_days = sum(fifelse(wet_day, get(value_col), 0), na.rm = TRUE),
      wet10mm_days = sum(wet10mm, na.rm = TRUE),
      wet20mm_days = sum(wet20mm, na.rm = TRUE)
    ),
    by = .(x, y, year)
  ]
  
  yearly <- by_year[
    ,
    .(
      avg_wet_days = mean(wet_days, na.rm = TRUE),
      avg_dry_days = mean(dry_days, na.rm = TRUE),
      avg_total_precip_wet_days = mean(total_precip_wet_days, na.rm = TRUE),
      avg_wet10mm_days = mean(wet10mm_days, na.rm = TRUE),
      avg_wet20mm_days = mean(wet20mm_days, na.rm = TRUE)
    ),
    by = .(x, y)
  ][, `:=`(freq = "yearly", period = 1)]
  
  
  # --- Combine ---
  combined <- rbindlist(list(monthly, quarterly, yearly), use.names = TRUE, fill = TRUE)
  
  # Join quantiles
  result <- merge(combined, quantiles, by = c("x", "y"), all.x = TRUE)
  
  return(list(
    result = result[]
    ))
  
}

#===============================================================================

# Spatial desagregation function at the communes level

spatial_agg_era5 <- function(
    source_dir,
    start_year,
    end_year,
    start_reference_year,
    end_reference_year
    ){
  
  source_dir <- dir$prepared
  start_year <- "2006"
  end_year <- "2023"
  start_reference_year <- "1971"
  end_reference_year <- "2000"
  
  # 1). Find all weather .rds file with right periods 
  files <- list.files(source_dir, full.names = TRUE, pattern = "\\.rds$")
  weather_files <- files[grepl("weather", files) & !(grepl("communes", files))]
  period_weather_files <- weather_files[grepl(paste0(start_year, "_", end_year), weather_files)]
  ref_period_weather_files <- period_weather_files[grepl(paste0(start_reference_year, "_", end_reference_year), period_weather_files)]
  
  # 2). Load ERA5 original raster file to extract CRS and have a raster template
  original_raster <- list.files(here(dir$source), full.names = TRUE, pattern = "\\.nc$")
  base_name <- tools::file_path_sans_ext(basename(original_raster))
  
  original_nc <- file.path(dir$source, paste0(base_name, ".nc"))
  rast_template <- rast(original_nc)
  ref_crs <- crs(rast_template)
  
  # 3). Get commune shapefile, transform to same CRS and get communes centroids
  communes <- st_read(here(dir$sf,"communes-20220101.shp"))
  # And convert to same crs as the template
  communes <- st_transform(communes, crs(ref_crs))
  # Get communes centroids
  communes_centroids <- communes %>% 
    st_centroid() %>%
    dplyr::select(insee, nom) 
  
  # 4). Convert the .rds file back to raster using the original file CRS
  for (file in ref_period_weather_files) {
    
    print(paste("Processing file:", file))
    
    # file <- ref_period_weather_files[1]
    df <- readRDS(file)
    
    # Weather unique cells
    weather_unique_cells <- df %>%
      distinct(x, y) %>%
      mutate(cell_id = row_number())
    
    # Get raster points from df
    raster_points <- weather_unique_cells %>%
      st_as_sf(coords = c("x", "y"), crs = ref_crs)
    
    # Extract coords as matrix for FNN
    raster_coords <- st_coordinates(raster_points)
    commune_coords <- st_coordinates(communes_centroids)
    
    # Find 4 nearest raster points from each commune centroids
    k <- 4  # how many neighbors
    nn <- get.knnx(raster_coords, commune_coords, k = k)
    
    # For each commune: get neighbor indices and distances
    neighbors_df <- data.frame(
      insee = rep(communes_centroids$insee, each = k),
      name       = rep(communes_centroids$nom, each = k),
      cell_id = as.vector(t(nn$nn.index)),
      distance =  as.vector(t(nn$nn.dist))
    )
    # Compute weights
    neighbors_df <- neighbors_df %>% 
      group_by(insee) %>% 
      mutate(
        weight = (1 / distance) / (sum(1 / distance))
      )
    
    # Join back (x,y) using cell_id
    neighbors_df <- neighbors_df %>%
      left_join(
        weather_unique_cells, 
        by = "cell_id"
      )
    
    # Join with raster data
    neighbors_with_data <- neighbors_df %>%
      left_join(
        df,
        by = c("x", "y"), 
        relationship = "many-to-many"
      ) %>% 
      ungroup()
    
    # Convert to data.table if not already
    setDT(neighbors_with_data)
    
    # Columns to EXCLUDE from aggregation
    non_weather_vars <- c(
      "x", "y", "cell_id", "distance", "weight", "reference_id", "period"
    )
    
    # Get weather column names (everything except grouping vars and excluded vars)
    grouping_vars <- c("insee", "name", "year", "month", "quarter", "freq", "variable", "statistic", "reference_period")
    weather_vars <- setdiff(names(neighbors_with_data), c(grouping_vars, non_weather_vars))
    
    # data.table aggregation - MUCH faster!
    result <- neighbors_with_data[, 
                                  lapply(.SD, function(z) sum(z * weight, na.rm = TRUE)),
                                  by = .(insee, name, year, month, quarter, freq, variable, statistic, reference_period),
                                  .SDcols = weather_vars
    ]
      
    # Save results in a separate file for monthly, quarterly and yearly frequency
    result_yearly <- result %>%
      filter(freq == "yearly") %>% 
      as_tibble()
    result_quarterly <- result %>% 
      filter(freq == "quarterly") %>% 
      as_tibble()
    result_monthly <- result %>% 
      filter(freq == "monthly")%>% 
      as_tibble()
    
    # Extract variable and statistic
    variable_name <- result %>% 
      as_tibble() %>% 
      dplyr::select(variable) %>% 
      slice(1) %>% 
      pull()
    statistic_name <- result %>%
      as_tibble() %>% 
      dplyr::select(statistic) %>% 
      slice(1) %>% 
      pull()
    
    # Save results
    try(saveRDS(result_yearly, 
            file = here(
              dir$final, 
              paste0(
                "era5_",
                variable_name,
                "_",
                statistic_name,
                "_communes_yearly_",
                start_year, "_",
                end_year,
                "_reference_",
                start_reference_year, "_",
                end_reference_year,
                ".rds")
            )
    ))
    try(saveRDS(result_quarterly, 
            file = here(
              dir$final, 
              paste0(
                "era5_",
                variable_name,
                "_",
                statistic_name,
                "_communes_quarterly_",
                start_year, "_",
                end_year,
                "_reference_",
                start_reference_year, "_",
                end_reference_year,
                ".rds")
            )
    ))
    try(saveRDS(result_monthly, 
            file = here(
              dir$final, 
              paste0(
                "era5_",
                variable_name,
                "_",
                statistic_name,
                "_communes_monthly_",
                start_year, "_",
                end_year,
                "_reference_",
                start_reference_year, "_",
                end_reference_year,
                ".rds")
            )
    ))
    
    print(paste("File:", file, " processed and saved."))
    
    }
  

  
    }
