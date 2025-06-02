#===============================================================================
# Description: Functions to download and aggregate ERA5-LAND data
#===============================================================================

# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ecmwfr, tidyverse, dplyr, terra, maps, here, ncdf4, raster, climate, devtools, lubridate,
               sf, sp, rnaturalearth, Matrix, readr, zoo, gstat)

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

#===============================================================================

# Download function over multiple months for source data

download_era5_ecmwfr <- function(variable,
                                 start_date,
                                 end_date,
                                 statistic,
                                 area = c(52, -5, 42, 9),
                                 frequency = "6_hourly",
                                 time_zone = "utc+02:00",
                                 dataset = "derived-era5-single-levels-daily-statistics",
                                 product_type = "reanalysis",
                                 user = "ecmwfr",
                                 key,
                                 output_file = NULL,
                                 output_dir) {
  
  # 1). Set API key
  wf_set_key(key = key, user = user)
  
  # 2). Download separating by month
  
  # Generate dates
  dates <- seq(ymd(start_date), ymd(end_date), by = "day")
  date_df <- tibble(date = dates) |>
    mutate(year = format(date, "%Y"),
           month = format(date, "%m"),
           day = format(date, "%d")) |>
    group_by(year, month)
  
  downloaded_files <- c()
  
  # Split month by month
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
    
    # 3). Construct API request
    request <- list(
      dataset_short_name = dataset,
      product_type = product_type,
      variable = variable,
      year = year,
      month = month,
      day = days,
      daily_statistic = statistic,
      time_zone = time_zone,
      frequency = frequency,
      area = area,
      format = "netcdf",
      target = output_file_month
    )
    
    # 4). Launch request
    file_path <- wf_request(
      request = request,
      transfer = TRUE,
      user = user, 
      path = output_dir
    )
    
    downloaded_files <- c(downloaded_files, file_path)
  }
  
  return(invisible(downloaded_files))
}

#===============================================================================


# Automation of downloaded data processing and border delimitation

process_era5_files <- function(source_dir) {

  # 1). Find all .nc file for a given month and load them
  nc_files <- list.files(source_dir, pattern = "\\.nc$", full.names = TRUE)
  
  # 2). Crop raster data to France borders
  france <- ne_countries(scale = "medium", country = "France", returnclass = "sf")
  
  # 3). Transform into data frame format for each file
  processed_list <- lapply(nc_files, function(nc_file) {

    # Load netCDF file
    rast_data <- rast(nc_file)
    
    # Harmonize projection
    france_proj <- st_transform(france, crs(rast_data))
    
    # Crop data using France borders
    cropped <- crop(rast_data, vect(france_proj))
    masked <- mask(cropped, vect(france_proj))
    
    # Transform into data frame format
    df <- as.data.frame(masked, xy = TRUE, na.rm = TRUE)
    
    # Add treated filename
    df$file <- basename(nc_file)
    
    return(df)
  })
  
  # 4). Return all dataframes
  return(processed_list)
  
}

#===============================================================================


# Creation of basic time indicators : mean, minimum, maximum 

weather_df_indice <- function(df, prefixe_variable = c("t2m", "tp"), output_dir) {
  prefixe <- match.arg(prefixe_variable)
  
  # 1). Rename columns from 1 to n (number of columns in data frame)
  noms_colonnes <- colnames(df)
  
  # Detect columns with the associated pattern
  colonnes_valide <- grep(paste0("^", prefixe, "_valid_time=\\d+$"), noms_colonnes, value = TRUE)
  
  # Extract corresponding column numbers and replace name with it
  for (old_name in colonnes_valide) {
    i <- as.numeric(sub(paste0(prefixe, "_valid_time="), "", old_name))
    new_name <- as.character(i + 1)
    noms_colonnes[noms_colonnes == old_name] <- new_name
  }
  colnames(df) <- noms_colonnes
  
  # 2). Identify day columns
  colonnes_jours <- intersect(as.character(1:31), colnames(df))  ### can also be used as a robustness check
  
  # 3). Compute indicators
  df <- df |>
    mutate(
      moyenne = rowMeans(across(all_of(colonnes_jours)), na.rm = TRUE),
      minimum = apply(across(all_of(colonnes_jours)), 1, min, na.rm = TRUE),
      maximum = apply(across(all_of(colonnes_jours)), 1, max, na.rm = TRUE)
    )
  
  # 4). Save file in data -> prepared directory
  
  original_filename <- unique(df$file)
  file_name <- sub("\\.nc$", ".rds", original_filename)
  saveRDS(df, file.path(output_dir, file_name))
  
  return(df)
}

#===============================================================================

# Warm spell duration index

warm_spells_indices <- function(source_dir = dir$prepared,
                                normals_file = "data/prepared/normals.rds") {
  
  # 1). Get season normals file with the percentiles
  normals <- readRDS(normals_file)
  
  # 2). Find all .rds file for a given month and load them 
  files <- list.files(source_dir, full.names = TRUE, pattern = "\\.rds$")
  files <- files[basename(files) != "normals.rds"]
  
  # 3). Create warm spell and warm spell duration index for each file
  for (file in files) {
    
    df <- readRDS(file) 
    
    # Detect temperature columns (from 1 to n day in a month)
    jour_cols <- names(df)[str_detect(names(df), "^\\d+$")]
    jours_num <- as.integer(jour_cols)
    
    # Associate each date with the day of a year (like in "normals.rds")
    date_str <- str_extract(basename(file), "\\d{4}-\\d{2}")
    date_base <- as.Date(paste0(date_str, "-01"))
    day_base <- yday(date_base)  
    day_seq <- day_base + 0:(length(jour_cols) - 1)  
    
    # Extract the percentile columns of the month of interest
    p90_cols <- paste0("p90_day_", day_seq)
    
    df_p90 <- normals |> 
      dplyr::select(x, y, all_of(p90_cols))
    names(df_p90)[3:ncol(df_p90)] <- jour_cols ### rename columns as in temperature file 
    
    # Join both datasets
    df_full <- left_join(df, df_p90, by = c("x", "y"), suffix = c("", "_p90"))
    
    df_full <- df_full |> 
      rowwise() |> 
      mutate(
      # Compute warm spell = binary indicator that tells if there was a warm spell in the month
      warm_spell = {
        temp <- as.numeric(c_across(all_of(jour_cols)))
        seuil <- as.numeric(c_across(paste0(jour_cols, "_p90")))
        chaud <- temp > seuil
        as.integer(any(rollapply(chaud, 6, all, fill = FALSE, align = "left")))
      },
      # Compute warm spell duration = number of days under the warm spell
      warm_spell_index = {
        temp <- as.numeric(c_across(all_of(jour_cols)))
        seuil <- as.numeric(c_across(paste0(jour_cols, "_p90")))
        chaud <- temp > seuil
        indic <- rep(FALSE, length(chaud))
        pos <- which(rollapply(chaud, 6, all, fill = FALSE, align = "left"))
        for (i in pos) indic[i:(i + 5)] <- TRUE
        sum(indic)
      }
    ) |> ungroup()
    
    # Remove percentile columns
    df_final <- df_full |> 
      dplyr::select(-ends_with("_p90"))
    
    saveRDS(df_final, file)
  }
}

#Cold Spell duration index

cold_spells_indices <- function(source_dir = dir$prepared,
                                normals_file = "data/prepared/normals.rds") {
  
  # 1). Get season normals file with the percentiles
  normals <- readRDS(normals_file)
  
  # 2). Find all .rds file for a given month and load them 
  files <- list.files(source_dir, full.names = TRUE, pattern = "\\.rds$")
  files <- files[basename(files) != "normals.rds"]
  
  # 3). Create warm spell and warm spell duration index for each file
  for (file in files) {
    
    df <- readRDS(file) 
    
    # Detect temperature columns (from 1 to n day in a month)
    jour_cols <- names(df)[str_detect(names(df), "^\\d+$")]
    jours_num <- as.integer(jour_cols)
    
    # Associate each date with the day of a year (like in "normals.rds")
    date_str <- str_extract(basename(file), "\\d{4}-\\d{2}")
    date_base <- as.Date(paste0(date_str, "-01"))
    day_base <- yday(date_base)  
    day_seq <- day_base + 0:(length(jour_cols) - 1)  
    
    # Extract the percentile columns of the month of interest
    p10_cols <- paste0("p10_day_", day_seq)
    
    df_p10 <- normals |> 
      dplyr::select(x, y, all_of(p10_cols))
    names(df_p10)[3:ncol(df_p10)] <- jour_cols ### rename columns as in temperature file 
    
    # Join both datasets
    df_full <- left_join(df, df_p10, by = c("x", "y"), suffix = c("", "_p10"))
    
    df_full <- df_full |>
      rowwise() |> 
      mutate(
        # Compute cold spell = binary indicator that tells if there was a cold spell in the month
        cold_spell = {
          temp <- as.numeric(c_across(all_of(jour_cols)))
          seuil <- as.numeric(c_across(paste0(jour_cols, "_p10")))
          froid <- temp > seuil
          as.integer(any(rollapply(froid, 6, all, fill = FALSE, align = "left")))
        },
        # Compute cold spell duration = number of days under the cold spell
        cold_spell_index = {
          temp <- as.numeric(c_across(all_of(jour_cols)))
          seuil <- as.numeric(c_across(paste0(jour_cols, "_p10")))
          froid <- temp > seuil
          indic <- rep(FALSE, length(froid))
          pos <- which(rollapply(froid, 6, all, fill = FALSE, align = "left"))
          for (i in pos) indic[i:(i + 5)] <- TRUE
          sum(indic)
        }
      ) |> ungroup()
    
    # Remove percentile columns
    df_final <- df_full |> 
      dplyr::select(-ends_with("_p10"))
    
    saveRDS(df_final, file)
  }
}

#===============================================================================

# Spatial desagregation function

spatial_agg_era5 <- function(source_dir){
  
  # 1). Find all .rds file for a given month and load them 
  files <- list.files(source_dir, full.names = TRUE, pattern = "\\.rds$")
  files <- files[basename(files) != "normals.rds"]
  
  # 2). Get commune shapefile
  communes <- st_read("data/source/communes-20220101.shp")
  
  for (file in files) {
    
    # 3). Find for each .rds file the corresponding original raster file
    base_name <- tools::file_path_sans_ext(basename(file))
    
    original_nc <- file.path(dir$source, paste0(base_name, ".nc"))
    original_rast <- rast(original_nc)
    
    df <- readRDS(file)
    
    # 4). Convert the .rds file back to raster using the original file CRS
    spat <- vect(df, geom = c("x", "y"), crs = crs(original_rast))
    
    rast_template <- terra::rast(original_rast) ## according to ERA5 documentation
    
    value_cols <- setdiff(names(df), c("x", "y"))  # remove coord columns if needed
    
    rast_list <- lapply(value_cols, function(col) {
      rasterize(spat, rast_template, field = col, fun = "mean")
    })
    
    rast <- rast(rast_list)
    names(rast) <- value_cols
    
    # 5). Create a spatial vector of each communes centroid point
    
    # Transform the commune shapefile into a sf object
    communes <- st_transform(communes, crs(rast))
    
    # Get the communes centroid point
    communes$centre <- st_centroid(communes$geometry)
    
    # Extract (x, y) coordinates from centroid point
    coords <- st_coordinates(communes$centre)
    
    # Add (x, y) coordinates to the sf object
    communes$x <- coords[, 1]
    communes$y <- coords[, 2]
    
    # Convert centroid to spatial vector
    vect_centroid <- vect(communes$centre)
    
    # 6). Create data frame of spatial aggregation
    
    # Extract centroid values from weather raster
    rast_centroid <- extract(rast, vect_centroid)
    
    # Add commune coordinates to the created data frame
    final_centroid <- cbind(communes[, c("nom", "x", "y")], rast_centroid)
    
    # 7). Interpolate missing values with the Inverse Weighted Distance method
    idw_interpolate <- function(col_name) {
      col_vals <- final_centroid[[col_name]]
      # Skip columns that do not have to be interpolated
      if (all(is.na(col_vals))) {
        return(col_vals)
      }
      
      df <- data.frame(x = final_centroid$x,
                       y = final_centroid$y,
                       value = col_vals)
      
      # (Re)transform to dataframe/numeric to avoid eventual bugging
      df <- as.data.frame(df)
      df$x <- as.numeric(df$x)
      df$y <- as.numeric(df$y)
      
      df_valid <- df[!is.na(df$value), ]
      df_valid <- as.data.frame(df_valid)
      df_valid$x <- as.numeric(df_valid$x)
      df_valid$y <- as.numeric(df_valid$y)
      
      # Convert to spatial objects
      coordinates(df_valid) <- ~x + y
      coordinates(df) <- ~x + y
      
      # IDW
      model <- gstat(formula = value ~ 1, data = df_valid, nmax = 10, set = list(idp = 2))
      prediction <- predict(model, newdata = df)
      
      return(prediction$var1.pred)
    }
    
    # Remove unused columns from the data frame 
    final_centroid <- final_centroid |>
      select(-ID, -file)
    
    # Extract columns that need interpolating
    cols_static <- c("x", "y", "nom", "geometry")
    cols_to_interp <- setdiff(names(final_centroid), cols_static)
    
    # Apply interpolation to all numeric columns of the dataframe
    final_centroid_interp <- final_centroid
    final_centroid_interp[cols_to_interp] <- lapply(cols_to_interp, idw_interpolate)
    
    # Save final prepared dataset
    saveRDS(final_centroid_interp, file)
    }
}

