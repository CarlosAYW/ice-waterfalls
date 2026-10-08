# =====================================================================
# Copyright / usage information
# ---------------------------------------------------------------------
# 2025 All rights reserved.
#
# This script was created for the master thesis on icefall modeling
# in North Tyrol. The code may not be copied, modified, shared, or
# embedded in other projects (especially commercial or proprietary
# software) without the author's prior written permission.
#
# Usage beyond reading or reproducing this thesis context is prohibited.
# Please contact the author before requesting reuse.
# =====================================================================

suppressPackageStartupMessages({
  library(httr)
  library(raster)
  library(leaflet)
  library(htmlwidgets)
  library(htmltools)
  library(ncdf4)
  library(readr)
  library(dplyr)
  library(tibble)
  library(lubridate)
  library(png)
})

# 0) Time range & area -------------------------------------------------

start_all    <- as.Date("2025-10-01") # Season start date
end_all      <- Sys.Date()           # through today
chunk_days   <- 2                    # 2-day chunks (API limit)
chunk_starts <- seq.Date(start_all, end_all, by = chunk_days)

# North Tyrol bounding box in WGS84
bbox <- c(
  46.7,  # lat_min
  10.1,  # lon_min
  47.7,  # lat_max
  12.2   # lon_max
)

# INCA parameters (extended for wind, humidity, etc.)
parameters    <- c("RR", "T2M", "RH2M", "UU", "VV", "GL", "P0", "TD2M")
param_str     <- paste(parameters, collapse = ",")
base_url_inca <- "https://dataset.api.hub.geosphere.at/v1/grid/historical/inca-v1-1h-1km"

out_dir <- "data/inca_nordtirol"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

nc_files <- character(length(chunk_starts))

# 1) Load INCA data in chunks (skip existing files) --------------------

for (idx in seq_along(chunk_starts)) {
  cs <- as.Date(chunk_starts[idx])
  ce <- as.Date(min(cs + (chunk_days - 1), end_all))
  
  message("Chunk: ", cs, " to ", ce)
  
  start_time <- sprintf("%sT00:00", format(cs, "%Y-%m-%d"))
  end_time   <- sprintf("%sT23:00", format(ce, "%Y-%m-%d"))
  
  outfile <- file.path(
    out_dir,
    sprintf("inca_nordtirol_%s_%s.nc",
            format(cs, "%Y%m%d"),
            format(ce, "%Y%m%d"))
  )
  
  if (file.exists(outfile) && file.info(outfile)$size > 0) {
    message("  -> skipped (exists): ", outfile)
  } else {
    query <- list(
      parameters    = param_str,
      start         = start_time,
      end           = end_time,
      bbox          = paste(bbox, collapse = ","),
      output_format = "netcdf",
      filename      = "inca_nordtirol"
    )
    
    resp <- GET(
      url   = base_url_inca,
      query = query,
      write_disk(outfile, overwrite = TRUE)
    )
    stop_for_status(resp)
    
    message("  -> saved: ", outfile, " (", file.info(outfile)$size, " Bytes)")
  }
  
  nc_files[idx] <- outfile
}

nc_files <- sort(unique(nc_files))

# 2) Read all NetCDFs and pack them into one object --------------------

convert_nc_time <- function(time_vals, time_units) {
  if (!grepl("since", time_units)) {
    origin <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
    return(origin + time_vals * 3600)
  }
  unit_str   <- sub(" since.*", "", time_units)
  origin_str <- sub(".*since ", "", time_units)
  origin     <- as.POSIXct(origin_str, tz = "UTC")
  
  if (grepl("hour", unit_str, ignore.case = TRUE)) {
    origin + time_vals * 3600
  } else if (grepl("second", unit_str, ignore.case = TRUE)) {
    origin + time_vals
  } else if (grepl("day", unit_str, ignore.case = TRUE)) {
    origin + time_vals * 86400
  } else {
    origin + time_vals
  }
}

vars_wanted <- c("T2M", "RH2M", "RR", "GL", "UU", "VV", "TD2M", "P0")

nc0       <- nc_open(nc_files[1])
varnames0 <- names(nc0$var)

lon_name <- intersect(c("lon", "longitude", "x", "xc"), varnames0)[1]
lat_name <- intersect(c("lat", "latitude", "y", "yc"), varnames0)[1]

lon0 <- ncvar_get(nc0, lon_name)  # [nx, ny]
lat0 <- ncvar_get(nc0, lat_name)  # [nx, ny]

nx <- dim(lon0)[1]
ny <- dim(lon0)[2]

time_dim_name <- names(nc0$dim)[grep("time", tolower(names(nc0$dim)))[1]]
vars_found    <- intersect(vars_wanted, varnames0)
nc_close(nc0)

# Time lengths per file
nt_per_file <- integer(length(nc_files))
for (i in seq_along(nc_files)) {
  nc <- nc_open(nc_files[i])
  nt_per_file[i] <- length(ncvar_get(nc, time_dim_name))
  nc_close(nc)
}
nt_total <- sum(nt_per_file)

# Allocate arrays

time_all  <- as.POSIXct(rep(NA_real_, nt_total), origin = "1970-01-01", tz = "UTC")
inca_data <- lapply(vars_found, function(.) array(NA_real_, dim = c(nx, ny, nt_total)))
names(inca_data) <- vars_found

# Read data
pos <- 1L
for (i in seq_along(nc_files)) {
  nc <- nc_open(nc_files[i])
  
  time_vals  <- ncvar_get(nc, time_dim_name)
  time_units <- ncatt_get(nc, time_dim_name, "units")$value
  nt_i       <- length(time_vals)
  
  idx <- pos:(pos + nt_i - 1L)
  time_all[idx] <- convert_nc_time(time_vals, time_units)
  
  for (v in vars_found) {
    inca_data[[v]][ , , idx] <- ncvar_get(nc, v)
  }
  
  nc_close(nc)
  pos <- pos + nt_i
}

inca_nordtirol_all <- list(
  lon  = lon0,
  lat  = lat0,
  time = time_all,
  data = inca_data
)

# 3) FDH/MDH (hourly) + orientation + raster template ------------------

T2M_arr <- inca_nordtirol_all$data$T2M   # °C [nx, ny, nt]
RH_arr  <- inca_nordtirol_all$data$RH2M  # %  [nx, ny, nt]
UU_arr  <- inca_nordtirol_all$data$UU    # m/s
VV_arr  <- inca_nordtirol_all$data$VV    # m/s

lon <- inca_nordtirol_all$lon
lat <- inca_nordtirol_all$lat

nx <- dim(T2M_arr)[1]
ny <- dim(T2M_arr)[2]
nt <- dim(T2M_arr)[3]

# Orientation (S->N, W->E)
lat_mean_j <- colMeans(lat, na.rm = TRUE)
if (lat_mean_j[1] < tail(lat_mean_j, 1)) {
  idx_j <- ny:1
  lon   <- lon[, idx_j]
  lat   <- lat[, idx_j]
  T2M_arr <- T2M_arr[, idx_j, , drop = FALSE]
  RH_arr  <- RH_arr[,  idx_j, , drop = FALSE]
  UU_arr  <- UU_arr[,  idx_j, , drop = FALSE]
  VV_arr  <- VV_arr[,  idx_j, , drop = FALSE]
}

lon_mean_i <- rowMeans(lon, na.rm = TRUE)
if (lon_mean_i[1] > tail(lon_mean_i, 1)) {
  idx_i <- nx:1
  lon   <- lon[idx_i, ]
  lat   <- lat[idx_i, ]
  T2M_arr <- T2M_arr[idx_i, , , drop = FALSE]
  RH_arr  <- RH_arr[idx_i, , , drop = FALSE]
  UU_arr  <- UU_arr[idx_i, , , drop = FALSE]
  VV_arr  <- VV_arr[idx_i, , , drop = FALSE]
}

nx <- dim(T2M_arr)[1]
ny <- dim(T2M_arr)[2]
nt <- dim(T2M_arr)[3]

FDH_hourly <- pmax(-T2M_arr, 0)
MDH_hourly <- pmax( T2M_arr, 0)
W_arr      <- sqrt(UU_arr^2 + VV_arr^2)

wind_ref <- mean(W_arr, na.rm = TRUE)
if (!is.finite(wind_ref) || wind_ref <= 0) wind_ref <- 5

FDH_mat_plain <- t(apply(FDH_hourly, c(1, 2), sum, na.rm = TRUE)) # [ny, nx]

r_template <- raster(FDH_mat_plain)
extent(r_template) <- c(min(lon, na.rm = TRUE),
                        max(lon, na.rm = TRUE),
                        min(lat, na.rm = TRUE),
                        max(lat, na.rm = TRUE))
crs(r_template) <- "EPSG:4326"

# 4) DEM on INCA grid + exposure index --------------------------------

dem_inca <- raster("data/DEM/DEM_Tirol_INCAgrid_1km_epsg4326.tif")
crs(dem_inca) <- "EPSG:4326"
names(dem_inca) <- "elev_m"

sl_as       <- terrain(dem_inca, opt = c("slope", "aspect"), unit = "degrees")
aspect_inca <- sl_as[["aspect"]]

aspect_rad_from_north <- (aspect_inca - 90) * pi / 180
northness_r           <- cos(aspect_rad_from_north)

solar_index_r  <- (1 - northness_r) / 2
solar_index_r[is.na(solar_index_r[])] <- 0.5
solar_index_ij <- t(as.matrix(solar_index_r))

# --- Elevation weights: cells > 2500 m are weighted less --------------
alt_threshold    <- 2500
alt_weight_high  <- 0.5

alt_weight_r <- raster::calc(dem_inca, fun = function(z) {
  ifelse(is.na(z), NA,
         ifelse(z > alt_threshold, alt_weight_high, 1))
})
names(alt_weight_r) <- "w_alt"

weighted_mean_alt <- function(r, w_rast = alt_weight_r) {
  v <- raster::getValues(r)
  w <- raster::getValues(w_rast)
  ok <- is.finite(v) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(v[ok] * w[ok]) / sum(w[ok])
}

# 5) Time-dependent solar elevation + time weighting -------------------

time_vec <- inca_nordtirol_all$time
if (length(time_vec) != nt) {
  if (length(time_vec) > nt) time_vec <- tail(time_vec, nt)
  else stop("Length of time_vec (", length(time_vec), ") != nt (", nt, ")")
}

time_local <- as.POSIXlt(time_vec, tz = "Europe/Vienna")

doy  <- time_local$yday + 1
hour <- time_local$hour + time_local$min / 60 + time_local$sec / 3600

delta_t <- 23.44 * pi/180 * sin(2 * pi * (284 + doy) / 365)

lat_center_deg <- (min(lat, na.rm = TRUE) + max(lat, na.rm = TRUE)) / 2
lat_center_rad <- lat_center_deg * pi / 180

H_t <- (hour - 12) * 15 * pi / 180

sin_alpha_t <- sin(lat_center_rad) * sin(delta_t) +
  cos(lat_center_rad) * cos(delta_t) * cos(H_t)

sin_alpha_t[sin_alpha_t < 0] <- 0
solar_height_factor_t <- sin_alpha_t

# For ice-thickness accumulation, all hours are weighted equally.
# Time decay would artificially devalue early-season phases.
weight_time_t <- rep(1, length(time_vec))

t_start           <- min(time_vec, na.rm = TRUE)
time_offset_hours <- as.numeric(difftime(time_vec, t_start, units = "hours"))

# Local time for historical snapshots

time_local_all <- as.POSIXct(format(time_vec, tz = "Europe/Vienna", usetz = TRUE),
                             tz = "Europe/Vienna")

t_min_local <- min(time_local_all, na.rm = TRUE)
t_max_local <- max(time_local_all, na.rm = TRUE)

t0_local <- as.POSIXct(paste0(format(t_min_local, "%Y-%m-%d"), " 07:00:00"),
                       tz = "Europe/Vienna")
if (t0_local < t_min_local) t0_local <- t0_local + 24 * 3600

t_last_local <- as.POSIXct(paste0(format(t_max_local, "%Y-%m-%d"), " 07:00:00"),
                           tz = "Europe/Vienna")
if (t_last_local > t_max_local) t_last_local <- t_last_local - 24 * 3600

if (t_last_local < t0_local) {
  snap_hours_hist <- numeric(0)
} else {
  snap_times_local <- seq(from = t0_local, to = t_last_local, by = "1 day")
  snap_times_utc   <- as.POSIXct(format(snap_times_local, tz = "UTC", usetz = TRUE),
                                 tz = "UTC")
  snap_hours_hist  <- as.numeric(difftime(snap_times_utc, t_start, units = "hours"))
}

FDH_hist_layers <- list()
FDH_hist_labels <- character(0)
FDH_hist_times  <- as.POSIXct(character(0), tz = "UTC")
snap_idx_hist   <- 1L

# 6) Effective FDH (FDHm + wind + humidity + radiation + time) ---------

k_expo   <- 0.5
k_wind   <- 0.5
wind_min <- 0.5
wind_max <- 2.0
rh_opt   <- 0.7
rh_sig   <- 0.15
k_melt   <- 1.2

FDH_sum_eff <- matrix(0, nrow = nx, ncol = ny)

for (k in seq_len(nt)) {
  T_k  <- T2M_arr[ , , k]
  RH_k <- RH_arr[ , , k] / 100
  
  if (all(is.na(T_k)) || all(is.na(RH_k))) {
    t_k_hours <- time_offset_hours[k]
    while (snap_idx_hist <= length(snap_hours_hist) &&
           t_k_hours >= snap_hours_hist[snap_idx_hist] - 0.5) {
      
      FDH_k_snap <- FDH_sum_eff
      FDH_k_snap[FDH_k_snap < 0] <- 0
      
      r_FDH_snap <- raster(t(FDH_k_snap))
      extent(r_FDH_snap) <- extent(r_template)
      crs(r_FDH_snap)    <- crs(r_template)
      
      label_time <- t_start + snap_hours_hist[snap_idx_hist] * 3600
      label_str  <- format(label_time, "%d.%m.%Y")
      
      FDH_hist_layers[[length(FDH_hist_layers) + 1L]] <- r_FDH_snap
      FDH_hist_labels <- c(FDH_hist_labels, label_str)
      FDH_hist_times  <- c(FDH_hist_times, label_time)
      
      snap_idx_hist <- snap_idx_hist + 1L
    }
    next
  }
  
  RH_k_eff <- RH_k
  RH_k_eff[is.na(RH_k_eff)] <- rh_opt
  
  FDH_k <- pmax(-T_k, 0)
  MDH_k <- pmax( T_k, 0)
  W_k   <- W_arr[ , , k]
  
  wind_norm <- (W_k - wind_ref) / wind_ref
  f_wind    <- 1 + k_wind * wind_norm
  f_wind[!is.finite(f_wind)] <- 1
  f_wind <- pmax(wind_min, pmin(wind_max, f_wind))
  
  f_rh_peak <- exp(- (RH_k_eff - rh_opt)^2 / (2 * rh_sig^2))
  f_rh      <- 0.5 + 0.5 * f_rh_peak
  
  s_height <- solar_height_factor_t[k]
  if (s_height == 0) {
    f_rad <- 1
  } else {
    f_rad <- 1 - k_expo * s_height * solar_index_ij
  }
  f_rad[!is.finite(f_rad)] <- 1
  f_rad <- pmax(0, pmin(1, f_rad))
  
  f_time <- weight_time_t[k]
  if (!is.finite(f_time)) f_time <- 1
  
  FDH_eff_k <- FDH_k * f_wind * f_rh * f_rad * f_time
  MDH_eff_k <- MDH_k * k_melt * f_time
  
  FDH_sum_eff <- FDH_sum_eff + (FDH_eff_k - MDH_eff_k)
  
  t_k_hours <- time_offset_hours[k]
  while (snap_idx_hist <= length(snap_hours_hist) &&
         t_k_hours >= snap_hours_hist[snap_idx_hist] - 0.5) {
    
    FDH_k_snap <- FDH_sum_eff
    FDH_k_snap[FDH_k_snap < 0] <- 0
    
    r_FDH_snap <- raster(t(FDH_k_snap))
    extent(r_FDH_snap) <- extent(r_template)
    crs(r_FDH_snap)    <- crs(r_template)
    
    label_time <- t_start + snap_hours_hist[snap_idx_hist] * 3600
    label_str  <- format(label_time, "%d.%m.%Y")
    
    FDH_hist_layers[[length(FDH_hist_layers) + 1L]] <- r_FDH_snap
    FDH_hist_labels <- c(FDH_hist_labels, label_str)
    FDH_hist_times  <- c(FDH_hist_times, label_time)
    
    snap_idx_hist <- snap_idx_hist + 1L
  }
}

FDH_sum_eff[FDH_sum_eff < 0] <- 0

r_FDH_eff <- raster(t(FDH_sum_eff))
extent(r_FDH_eff) <- extent(r_template)
crs(r_FDH_eff)    <- crs(r_template)
names(r_FDH_eff)  <- "FDH_eff_C_h"
# 7) Effective FDH -> ice thickness -----------------------------------

h_c      <- 30
rho_i    <- 880
Lf       <- 334000
FDH_crit <- 1

alpha <- h_c * 3600 / (rho_i * Lf)

ice_hist_layers <- list()
ice_hist_labels <- character(0)

if (length(FDH_hist_layers) > 0) {
  ice_hist_layers <- vector("list", length(FDH_hist_layers))
  ice_hist_labels <- FDH_hist_labels
  for (i in seq_along(FDH_hist_layers)) {
    r_fd  <- FDH_hist_layers[[i]]
    ice_i <- r_fd * alpha
    ice_i[r_fd < FDH_crit] <- NA
    names(ice_i) <- paste0("h_pot_expo_hist_", ice_hist_labels[i])
    ice_hist_layers[[i]] <- ice_i
  }
}

ice_thick_expo <- r_FDH_eff * alpha
ice_thick_expo[r_FDH_eff < FDH_crit] <- NA
names(ice_thick_expo) <- "h_pot_expo_m"


# 9) Controls / HTML summary ------------------------------------------
# 10) Build time layers -----------------------------------------------

ice_time_layers <- list()
time_labels     <- character(0)

if (length(ice_hist_layers) > 0L) {
  ice_time_layers <- c(ice_time_layers, ice_hist_layers)
  time_labels     <- c(time_labels, ice_hist_labels)
}

if (length(ice_time_layers) == 0L) {
  ice_time_layers <- list(ice_thick_expo)
  time_labels     <- format(Sys.Date(), "%d.%m.%Y (Jetzt)")
}

n_steps <- min(length(ice_time_layers), length(time_labels))
ice_time_layers <- ice_time_layers[seq_len(n_steps)]
time_labels     <- time_labels[seq_len(n_steps)]

max_h_all <- vapply(
  ice_time_layers,
  function(r) {
    if (is.null(r)) return(NA_real_)
    suppressWarnings(cellStats(r, "max", na.rm = TRUE))
  },
  numeric(1)
)
max_h <- max(max_h_all, na.rm = TRUE)
if (!is.finite(max_h) || max_h <= 0) max_h <- 1

col_fun <- colorRampPalette(c("#ffffff", "#c6dbef", "#6baed6", "#08519c"))

pal_h <- colorNumeric(
  palette  = col_fun(100),
  domain   = c(0, max_h),
  na.color = "transparent"
)

last_update <- format(Sys.time(), tz = "Europe/Vienna", "%d.%m.%Y %H:%M %Z")
last_update_iso <- format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ")
last_update_payload <- jsonlite::toJSON(
  list(
    last_update = last_update,
    last_update_iso = last_update_iso
  ),
  auto_unbox = TRUE,
  pretty = TRUE
)
dir.create("site", showWarnings = FALSE)
writeLines(last_update_payload, "site/last_update.json", useBytes = TRUE)
writeLines(last_update_payload, "last_update.json", useBytes = TRUE)
ext         <- extent(r_template)

# =====================================================================
# Write external PNGs per step
# =====================================================================

dir.create("site/img", recursive = TRUE, showWarnings = FALSE)
dir.create("site/plots", recursive = TRUE, showWarnings = FALSE)

hex_to_rgba <- function(hex) {
  hex <- gsub("#", "", hex)
  hex[hex == "transparent" | is.na(hex) | nchar(hex) == 0] <- "00000000"
  hex[nchar(hex) == 6] <- paste0(hex[nchar(hex) == 6], "FF")
  r <- strtoi(substr(hex, 1, 2), 16L) / 255
  g <- strtoi(substr(hex, 3, 4), 16L) / 255
  b <- strtoi(substr(hex, 5, 6), 16L) / 255
  a <- strtoi(substr(hex, 7, 8), 16L) / 255
  cbind(r, g, b, a)
}

write_overlay_png <- function(r, pal, file) {
  v <- raster::getValues(r)
  cols <- pal(v)
  cols[!is.finite(v)] <- "#00000000"
  cols[cols == "transparent"] <- "#00000000"
  
  rgba <- hex_to_rgba(cols)
  nr <- raster::nrow(r)
  nc <- raster::ncol(r)
  
  arr <- array(0, dim = c(nr, nc, 4))
  arr[, , 1] <- matrix(rgba[,1], nr, nc, byrow = TRUE)
  arr[, , 2] <- matrix(rgba[,2], nr, nc, byrow = TRUE)
  arr[, , 3] <- matrix(rgba[,3], nr, nc, byrow = TRUE)
  arr[, , 4] <- matrix(rgba[,4], nr, nc, byrow = TRUE)
  
  png::writePNG(arr, target = file)
}

message("Writing PNG overlays for ", n_steps, " steps ...")
for (i in seq_len(n_steps)) {
  write_overlay_png(ice_time_layers[[i]], pal_h, sprintf("site/img/ice_%03d.png", i))
}

# 11) Icefall sun data + topo URLs ------------------------------------

DIR_SUN <- "data/suntime"

find_sun_file_for_uid <- function(uid, dir_sun = DIR_SUN) {
  uid_i <- suppressWarnings(as.integer(uid))
  if (!is.finite(uid_i)) return(NA_character_)
  cand <- c(
    file.path(dir_sun, sprintf("sun_uid_%03d.csv", uid_i)),
    file.path(dir_sun, sprintf("sun_uid_%d.csv", uid_i))
  )
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

load_sun_for_uids <- function(uids, dir_sun = DIR_SUN) {
  out <- vector("list", length(uids))
  missing_uids <- integer(0)
  for (i in seq_along(uids)) {
    uid <- as.integer(uids[[i]])
    f <- find_sun_file_for_uid(uid, dir_sun = dir_sun)
    if (is.na(f)) {
      missing_uids <- c(missing_uids, uid)
      next
    }
    df <- tryCatch(readr::read_csv(f, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) {
      missing_uids <- c(missing_uids, uid)
      next
    }
    out[[i]] <- df
  }
  list(
    data = dplyr::bind_rows(out),
    missing_uids = sort(unique(missing_uids))
  )
}

# Optional meta (difficulty) for map filters
PATH_META <- "data/Koordinaten_Wasserfaelle/eisklettern_links_entries_diff.csv"
meta_map <- NULL
if (file.exists(PATH_META)) {
  parse_uid <- function(x) {
    as.integer(readr::parse_number(as.character(x)))
  }
  get_chr <- function(df, ...) {
    cands <- c(...)
    for (nm in cands) if (nm %in% names(df)) return(as.character(df[[nm]]))
    rep(NA_character_, nrow(df))
  }
  to_num <- function(x) {
    if (is.null(x)) return(NA_real_)
    x <- as.character(x)
    x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
    x <- gsub(",", ".", x, fixed = TRUE)
    suppressWarnings(as.numeric(x))
  }
  meta_raw <- readr::read_delim(
    PATH_META,
    delim = ";",
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  ) %>%
    rename_with(tolower)
  meta_map <- tibble(
    uid = parse_uid(meta_raw$uid),
    name = get_chr(meta_raw, "name", "icefall", "icefall_name"),
    difficulty = get_chr(meta_raw, "difficulty", "grade"),
    elev_m = to_num(get_chr(meta_raw, "elevation_dgm5m", "elevation_m", "elevation", "elev_m", "altitude_m")),
    icefall_height_m = to_num(get_chr(meta_raw, "icefall_height_m", "height_m", "height")),
    topo_url = get_chr(meta_raw, "topo_url"),
    latitude = to_num(get_chr(meta_raw, "latitude")),
    longitude = to_num(get_chr(meta_raw, "longitude"))
  )
}

if (!dir.exists(DIR_SUN)) {
  message("Warning: sun directory missing: ", DIR_SUN)
}

expected_uids <- if (!is.null(meta_map) && "uid" %in% names(meta_map)) {
  sort(unique(meta_map$uid[is.finite(meta_map$uid)]))
} else {
  integer(0)
}

if (length(expected_uids) > 0) {
  if (dir.exists(DIR_SUN)) {
    sun_loaded <- load_sun_for_uids(expected_uids, dir_sun = DIR_SUN)
    if (length(sun_loaded$missing_uids) > 0) {
      message(
        "Warning: missing sun files for UIDs: ",
        paste(sprintf("%03d", sun_loaded$missing_uids), collapse = ", ")
      )
    }
    sun_df <- sun_loaded$data
  } else {
    message(
      "Warning: missing sun files for UIDs: ",
      paste(sprintf("%03d", expected_uids), collapse = ", ")
    )
    sun_df <- tibble::tibble()
  }
} else {
  if (dir.exists(DIR_SUN)) {
    sun_files <- sort(list.files(DIR_SUN, pattern = "^sun_uid_\\d+\\.csv$", full.names = TRUE))
    if (length(sun_files) == 0) {
      message("Warning: no sun files found in: ", DIR_SUN)
      sun_df <- tibble::tibble()
    } else {
      sun_df <- dplyr::bind_rows(lapply(sun_files, function(f) readr::read_csv(f, show_col_types = FALSE)))
    }
  } else {
    sun_df <- tibble::tibble()
  }
}

sun_missing_uids <- if (length(expected_uids) > 0) {
  setdiff(expected_uids, unique(sun_df$uid[is.finite(sun_df$uid)]))
} else {
  integer(0)
}

if (nrow(sun_df) == 0) {
  sun_df <- tibble::tibble(
    uid = integer(),
    date = as.Date(character()),
    sunrise_topo = as.POSIXct(character(), tz = "UTC"),
    sunset_topo = as.POSIXct(character(), tz = "UTC"),
    sun_hours_topo = numeric()
  )
}

sun_df <- sun_df %>%
  dplyr::mutate(
    uid  = as.integer(readr::parse_number(as.character(uid))),
    date = as.Date(date)
  )

if (!inherits(sun_df$sunrise_topo, "POSIXt")) {
  sun_df <- sun_df %>%
    mutate(
      sunrise_topo = ymd_hms(sunrise_topo, tz = "UTC"),
      sunset_topo  = ymd_hms(sunset_topo,  tz = "UTC")
    )
}

sun_df <- sun_df %>%
  mutate(
    sunrise_topo = with_tz(sunrise_topo, "Europe/Vienna"),
    sunset_topo  = with_tz(sunset_topo,  "Europe/Vienna"),
    sun_hours_topo = as.numeric(difftime(sunset_topo, sunrise_topo, units = "hours"))
  )

sun_date <- as.Date(Sys.time(), tz = "Europe/Vienna")

sun_today <- sun_df %>% dplyr::filter(date == sun_date)

if (nrow(sun_today) == 0) {
  last_date <- max(sun_df$date, na.rm = TRUE)
  message("No sun data found for ", sun_date, ". Using ", last_date, " instead.")
  sun_today <- sun_df %>% dplyr::filter(date == last_date)
}

if (nrow(sun_today) == 0) {
  message("Warning: sun_today is empty; no sun times will be shown.")
}

if (!is.null(meta_map)) {
  marker_base <- meta_map %>%
    dplyr::filter(is.finite(uid)) %>%
    dplyr::select(uid, name, latitude, longitude, topo_url, difficulty, elev_m, icefall_height_m) %>%
    dplyr::arrange(uid) %>%
    dplyr::distinct(uid, .keep_all = TRUE)

  if (nrow(marker_base) == 0) {
    message("Warning: metadata contains no valid UIDs; using sun data as marker base.")
    sun_today <- sun_today
  } else {
    sun_today <- marker_base %>%
      dplyr::left_join(
        sun_today %>%
          dplyr::select(uid, name, date, sunrise_topo, sunset_topo, sun_hours_topo),
        by = "uid",
        suffix = c("_meta", "_sun")
      ) %>%
      dplyr::mutate(
        name = dplyr::coalesce(.data$name_meta, .data$name_sun),
        date = dplyr::coalesce(.data$date, sun_date),
        sun_status = dplyr::case_when(
          uid %in% sun_missing_uids ~ "missing_file",
          is.finite(sun_hours_topo) & !is.na(sunrise_topo) & !is.na(sunset_topo) ~ "has_values",
          TRUE ~ "no_values_day"
        )
      ) %>%
      dplyr::select(-dplyr::any_of(c("name_meta", "name_sun")))
  }
} else {
  sun_today <- sun_today %>%
    dplyr::mutate(
      difficulty = NA_character_,
      elev_m = NA_real_,
      icefall_height_m = NA_real_,
      sun_status = dplyr::if_else(
        is.finite(sun_hours_topo) & !is.na(sunrise_topo) & !is.na(sunset_topo),
        "has_values",
        "no_values_day"
      )
    )
}

if (!"topo_url" %in% names(sun_today)) {
  sun_today$topo_url <- NA_character_
}
if (!"latitude" %in% names(sun_today)) {
  sun_today$latitude <- NA_real_
}
if (!"longitude" %in% names(sun_today)) {
  sun_today$longitude <- NA_real_
}

sun_today <- sun_today %>%
  dplyr::mutate(
    name = dplyr::coalesce(name, paste0("Icefall ", sprintf("%03d", uid))),
    sunrise_txt   = substr(as.character(sunrise_topo), 12, 16),
    sunset_txt    = substr(as.character(sunset_topo),  12, 16),
    sun_hours_txt = sprintf("%.1f", sun_hours_topo),
    date_txt      = format(date, "%d.%m.%Y"),
    link_txt = ifelse(
      !is.na(topo_url) & topo_url != "",
      paste0("<a href='", topo_url, "' target='_blank'>Open topo</a>"),
      "(no topo link available)"
    ),
    
    uid_pad  = sprintf("%03d", uid),
    plot_png = paste0("plots/uid_", uid_pad, ".png"),
    
    detail_url = paste0("icefalls/uid_", uid_pad, ".html"),

    difficulty_txt = dplyr::if_else(
      !is.na(difficulty) & difficulty != "",
      difficulty,
      "n/a"
    ),
    height_txt = dplyr::case_when(
      is.finite(icefall_height_m) ~ paste0(round(icefall_height_m), " m"),
      is.finite(elev_m) ~ paste0(round(elev_m), " m"),
      TRUE ~ "n/a"
    ),
    info_block = sprintf(
      "<div style='margin-top:4px;font-size:12px;line-height:1.35;'><b>Difficulty:</b> %s<br/><b>Elevation:</b> %s</div>",
      htmltools::htmlEscape(difficulty_txt),
      htmltools::htmlEscape(height_txt)
    ),

    map_meta = paste0(
      "<span class='map-meta' data-uid='", uid,
      "' data-name='", htmltools::htmlEscape(ifelse(is.na(name), "", name), attribute = TRUE),
      "' data-difficulty='", htmltools::htmlEscape(ifelse(is.na(difficulty), "", difficulty), attribute = TRUE),
      "' data-sun='", ifelse(is.na(sun_hours_topo), "", sprintf("%.2f", sun_hours_topo)),
      "'></span>"
    ),
    
    plot_block = paste0(
      "<hr style='margin:6px 0;'/>",
      
      "<div style='display:flex; gap:6px; flex-wrap:wrap;'>",
      
      "<a href='", detail_url, "' ",
      "style='padding:6px 10px; background:#0d6efd; color:white; ",
      "border-radius:6px; text-decoration:none; font-weight:600;'>",
      "📄 Details & upload",
      "</a>",

      "</div>",
      
      "<a href='", plot_png, "' target='_blank'>",
      "<img src='", plot_png, "' ",
      "style='width:320px;max-width:100%;height:auto;",
      "border:1px solid #ccc;border-radius:4px;margin-top:6px;' ",
      "onerror=\"this.style.display='none';\"/>",
      "</a>"
    ),
    
    popup = ifelse(
      sun_status == "missing_file",
      paste0(
        map_meta,
        sprintf(
          "<b>%s</b><br/>Sun on %s: NO SUN DATA%s<br/>%s",
          name, date_txt, info_block, link_txt
        ),
        plot_block
      ),
      ifelse(
        is.na(sun_hours_topo) | is.na(sunrise_topo) | is.na(sunset_topo),
        paste0(
          map_meta,
          sprintf(
          "<b>%s</b><br/>Sun on %s: no direct sunlight%s<br/>%s",
          name, date_txt, info_block, link_txt
          ),
          plot_block
        ),
        paste0(
          map_meta,
          sprintf(
            "<b>%s</b><br/>Sun on %s: %s – %s (%s h)%s<br/>%s",
            name, date_txt, sunrise_txt, sunset_txt, sun_hours_txt, info_block, link_txt
          ),
          plot_block
        )
      )
    )
  )

marker_data <- sun_today %>%
  dplyr::filter(is.finite(latitude), is.finite(longitude))

if (nrow(marker_data) == 0) {
  message("Warning: no markers with valid coordinates available; map shows no icefalls.")
} else {
  message("Markers ready: ", nrow(marker_data), " icefalls with valid coordinates.")
}

# 12) Leaflet map: stable ClusterGroup + slider swaps PNG URL ----------

init_i <- n_steps
m <- leaflet() |>
  addProviderTiles(providers$OpenStreetMap, group = "OSM") |>
  addProviderTiles(providers$OpenTopoMap,   group = "Terrain (Topo)")

preview_mode <- interactive()

if (isTRUE(preview_mode)) {
  m <- m |>
    addRasterImage(
      ice_time_layers[[init_i]],
      colors  = pal_h,
      opacity = 0.8,
      project = TRUE,
      method  = "bilinear",
      group   = "Ice thickness",
      layerId = "ice_preview"
    )
} else {
  bounds_js <- sprintf("[[%f,%f],[%f,%f]]", ext@ymin, ext@xmin, ext@ymax, ext@xmax)
  
  m <- htmlwidgets::onRender(
    m,
    sprintf(
      "function(el, x) {
         var map = this;
         var bounds = %s;
         function pad3(n){ return String(n).padStart(3,'0'); }

         var iceUrl = 'img/ice_' + pad3(%d) + '.png';

         var rasterPane = map.getPane('iceRasterPane');
         if (!rasterPane && typeof map.createPane === 'function') {
           rasterPane = map.createPane('iceRasterPane');
         }
         if (rasterPane && rasterPane.style) {
           rasterPane.style.zIndex = 250;
         }

         var ice = L.imageOverlay(iceUrl, bounds, {opacity: 0.8, layerId: 'ice_overlay', pane: 'iceRasterPane'});

         try {
           if (map.layerManager && typeof map.layerManager.addLayer === 'function') {
             map.layerManager.addLayer(ice, 'image', 'ice_overlay', 'Ice thickness', null, null);

             if (typeof map.layerManager.showGroup === 'function') {
               map.layerManager.showGroup('Ice thickness');
             }
           } else {
             ice.addTo(map);
           }
         } catch(e) {
           try { ice.addTo(map); } catch(e2) {}
         }

         window._iceOverlay = ice;
       }",
      bounds_js,
      init_i
    )
  )
}

# UI / controls
m <- m |>
  addMapPane("icefallsPane", zIndex = 650) |>
  addControl(
    position = "topleft",
    html = htmltools::HTML(
      "<div style='background:rgba(255,255,255,0.9);padding:6px 8px;border-radius:6px;display:flex;flex-direction:column;gap:6px;'>
         <a href='index.html' class='map-home-link' style='font-size:14px;font-weight:bold;'>🏠 Home</a>
         <a href='list.html' class='map-list-link' style='font-size:14px;font-weight:bold;'>📋 Icefall list</a>
       </div>"
    )
  ) |>
  addControl(
    position = "bottomright",
    html = htmltools::HTML(
      "<details id='map-filter' style='background:rgba(255,255,255,0.95);padding:8px 10px;border-radius:8px;box-shadow:0 4px 12px rgba(0,0,0,0.12);min-width:240px;width:340px;max-width:calc(100vw - 24px);box-sizing:border-box;'>
        <summary style='font-weight:700;font-size:12px;letter-spacing:0.02em;text-transform:uppercase;cursor:pointer;'>Filters</summary>
        <div style='margin-top:6px;display:flex;flex-direction:column;gap:10px;'>
          <input id='mapFilterInput' type='search' placeholder='Name, UID, difficulty' autocomplete='off'
                 style='width:100%;padding:6px 8px;border:1px solid #ddd;border-radius:8px;font-size:13px;'/>
          <div style='display:flex;flex-direction:column;gap:6px;'>
            <div style='font-size:11px;color:#444;font-weight:700;letter-spacing:0.02em;text-transform:uppercase;'>Difficulty</div>
            <div style='display:flex;flex-direction:column;gap:6px;'>
              <div style='display:flex;flex-direction:column;gap:4px;font-size:11px;color:#555;'>
                <span>Technical (A)</span>
                <div style='display:flex;gap:6px;align-items:center;flex-wrap:wrap;'>
                  <input id='mapAmin' type='range' min='0.75' max='4.25' step='0.25' value='0.75' style='flex:1;min-width:120px;'/>
                  <input id='mapAmax' type='range' min='0.75' max='4.25' step='0.25' value='4.25' style='flex:1;min-width:120px;'/>
                  <span id='mapARangeTxt' style='font-size:11px;color:#666;min-width:90px;'>A1- – A4+</span>
                </div>
              </div>

              <div style='display:flex;flex-direction:column;gap:4px;font-size:11px;color:#555;'>
                <span>Mixed (M)</span>
                <div style='display:flex;gap:6px;align-items:center;flex-wrap:wrap;'>
                  <input id='mapMmin' type='range' min='0.75' max='13.25' step='0.25' value='0.75' style='flex:1;min-width:120px;'/>
                  <input id='mapMmax' type='range' min='0.75' max='13.25' step='0.25' value='13.25' style='flex:1;min-width:120px;'/>
                  <span id='mapMRangeTxt' style='font-size:11px;color:#666;min-width:90px;'>M1- – M13+</span>
                </div>
              </div>

              <div style='display:flex;flex-direction:column;gap:4px;font-size:11px;color:#555;'>
                <span>Water ice (WI)</span>
                <div style='display:flex;gap:6px;align-items:center;flex-wrap:wrap;'>
                  <input id='mapWImin' type='range' min='0.75' max='7.25' step='0.25' value='0.75' style='flex:1;min-width:120px;'/>
                  <input id='mapWImax' type='range' min='0.75' max='7.25' step='0.25' value='7.25' style='flex:1;min-width:120px;'/>
                  <span id='mapWIRangeTxt' style='font-size:11px;color:#666;min-width:90px;'>WI1- – WI7+</span>
                </div>
              </div>

              <div style='display:flex;flex-direction:column;gap:4px;font-size:11px;color:#555;'>
                <span>Rock (UIAA)</span>
                <div style='display:flex;gap:6px;align-items:center;flex-wrap:wrap;'>
                  <input id='mapRmin' type='range' min='0.75' max='12.25' step='0.25' value='0.75' style='flex:1;min-width:120px;'/>
                  <input id='mapRmax' type='range' min='0.75' max='12.25' step='0.25' value='12.25' style='flex:1;min-width:120px;'/>
                  <span id='mapRRangeTxt' style='font-size:11px;color:#666;min-width:90px;'>1- – 12+</span>
                </div>
              </div>
            </div>
          </div>

          <div style='display:flex;flex-direction:column;gap:6px;'>
            <div style='font-size:11px;color:#444;font-weight:700;letter-spacing:0.02em;text-transform:uppercase;'>Sun</div>
            <div style='display:flex;flex-direction:column;gap:4px;font-size:11px;color:#555;'>
              <span>Sun today (h)</span>
              <div style='display:flex;gap:6px;align-items:center;flex-wrap:wrap;'>
                <input id='mapSunMin' type='range' min='0' max='12' step='0.25' value='0' style='flex:1;min-width:120px;'/>
                <input id='mapSunMax' type='range' min='0' max='12' step='0.25' value='12' style='flex:1;min-width:120px;'/>
                <span id='mapSunRangeTxt' style='font-size:11px;color:#666;min-width:90px;'>0.0 – 12.0 h</span>
              </div>
            </div>
          </div>

          <button id='mapFilterReset' type='button'
                  style='border:1px solid #d1d5db;background:#fff;border-radius:8px;padding:6px 8px;font-size:12px;cursor:pointer;'>
            Reset filters
          </button>
          <div id='mapFilterStatus' style='font-size:12px;color:#666;line-height:1.2;min-height:1.2em;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;width:100%;'></div>
        </div>
      </details>"
    )
  ) |>
  addControl(
    position = "bottomright",
    html = htmltools::HTML(
      "<button id='mapUseGeo' type='button' title='GPS'
              style='width:34px;height:34px;border-radius:999px;border:1px solid #d1d5db;background:#fff;box-shadow:0 4px 10px rgba(0,0,0,0.12);font-size:16px;cursor:pointer;margin-top:8px;'>
         📍
       </button>"
    )
  )

# Raw markers in their own group (not clustered in R)
m <- m |>
  addMarkers(
    data   = marker_data,
    lng    = ~longitude,
    lat    = ~latitude,
    popup  = ~popup,
    group  = "IcefallsRaw",
    options = markerOptions(pane = "icefallsPane"),
    layerId = ~uid
  )

# Dependency patch: load leaflet.markercluster even though R does not cluster the markers.
m <- m |>
  addMarkers(
    data = marker_data[0, , drop = FALSE],
    lng = ~longitude, lat = ~latitude,
    clusterOptions = markerClusterOptions()
  )

# Layers control (shows "Icefalls" as overlay and hides "IcefallsRaw")
m <- m |>
  addLayersControl(
    baseGroups    = c("OSM", "Terrain (Topo)"),
    overlayGroups = c("Ice thickness", "Icefalls"),
    options       = layersControlOptions(collapsed = FALSE)
  ) |>
  fitBounds(lng1 = ext@xmin, lat1 = ext@ymin, lng2 = ext@xmax, lat2 = ext@ymax)

# Legends
m <- m |>
  addLegend(
    pal       = pal_h,
    values    = c(0, max_h),
    title     = "Ice thickness (m)",
    labFormat = labelFormat(digits = 2),
    position  = "bottomleft"
  )

# --- Cluster + filter logic from scripts/map_cluster_filter.js (only clear/add on cluster)
js_cluster_filter <- paste(readLines("scripts/map_cluster_filter.js", warn = FALSE), collapse = "\n")

m <- htmlwidgets::onRender(m, js_cluster_filter)

# Time slider: steps stay 1:1, but only the PNG URL changes -----------
m <- m |>
  (\(x) {
    if (length(time_labels) > 0L) {
      labels_js  <- paste0("['", paste(time_labels, collapse = "','"), "']")
      n_steps_js <- n_steps

      js_code <- sprintf(
        "function(el, x) {
          var map = this;
          
          var isMobile = window.matchMedia && window.matchMedia('(max-width: 720px)').matches;
          
          var lc = el.getElementsByClassName('leaflet-control-layers-expanded')[0]
          || el.getElementsByClassName('leaflet-control-layers')[0];
          if (lc) {
            if (isMobile) {
              lc.style.marginTop   = '6px';
              lc.style.marginRight = '0';
              lc.style.transform   = 'none';
              lc.style.padding     = '8px 10px';
              lc.style.fontSize    = '13px';
            } else {
              lc.style.marginTop   = '10px';
              lc.style.marginRight = '90px';
              lc.style.transform   = 'scale(1.5)';
              lc.style.transformOrigin = 'top left';
              lc.style.padding     = '12px 15px';
              lc.style.fontSize    = '16px';
            }
          }
          
          var labels = %s;
          var nSteps = %d;
          if (!labels || labels.length === 0 || nSteps <= 0) return;
          
          function pad3(n){ return String(n).padStart(3,'0'); }
          
          var iceLayer = window._iceOverlay || null;
          
          if (!iceLayer) {
            map.eachLayer(function(l){
              if(!iceLayer && l && l.options && l.options.layerId === 'ice_overlay') iceLayer = l;
            });
          }
          
          var initial = nSteps - 1;
          if (initial < 0) initial = 0;
          
          function setTimeStep(step) {
            if (step < 0) step = 0;
            if (step >= nSteps) step = nSteps - 1;
            
            var i = step + 1; // 1..nSteps
            if (iceLayer) iceLayer.setUrl('img/ice_' + pad3(i) + '.png');
            
            var labelDiv = document.getElementById('time-label');
            if (labelDiv && step >= 0 && step < labels.length) {
              labelDiv.textContent = labels[step];
            }
            
            var next = Math.min(nSteps, i+1);
            var prev = Math.max(1, i-1);
            var img1 = new Image(); img1.src = 'img/ice_' + pad3(next) + '.png';
            var img3 = new Image(); img3.src = 'img/ice_' + pad3(prev) + '.png';
          }
          
          var sliderControl = L.control({position: 'topright'});
          sliderControl.onAdd = function() {
            var div = L.DomUtil.create('div', 'info leaflet-control');
            div.style.background   = 'rgba(255,255,255,0.9)';
            div.style.padding      = '8px 10px';
            div.style.borderRadius = '6px';
            div.style.minWidth     = isMobile ? '200px' : '260px';
            
            div.style.marginTop    = isMobile ? '110px' : '140px';
            div.style.marginRight  = isMobile ? '6px' : '10px';
            
            var title = document.createElement('div');
            title.style.fontSize    = '16px';
            title.style.marginBottom = '4px';
            title.innerHTML = '<b>Ice thickness timeline</b>';
            div.appendChild(title);
            
            var labelDiv = document.createElement('div');
            labelDiv.id = 'time-label';
            labelDiv.style.fontSize   = '14px';
            labelDiv.style.marginBottom = '4px';
            labelDiv.textContent = labels[initial] || labels[0];
            div.appendChild(labelDiv);
            
            var slider = document.createElement('input');
            slider.type  = 'range';
            slider.min   = 0;
            slider.max   = nSteps - 1;
            slider.step  = 1;
            slider.value = initial;
            slider.style.width = isMobile ? '200px' : '240px';
            slider.id    = 'time-slider';
            div.appendChild(slider);
            
            slider.addEventListener('input', function(e) {
              var step = parseInt(e.target.value, 10);
              if (!isNaN(step)) setTimeStep(step);
            });
            
            slider.addEventListener('mousedown', function() { if (map && map.dragging) map.dragging.disable(); });
            slider.addEventListener('mouseup',   function() { if (map && map.dragging) map.dragging.enable();  });
            
            L.DomEvent.disableClickPropagation(div);
            return div;
          };
          sliderControl.addTo(map);
          
          setTimeout(function(){ setTimeStep(initial); }, 200);
        }",
        labels_js,
        n_steps_js
      )

      htmlwidgets::onRender(x, js_code)
    } else {
      x
    }
  })()

# 13) Output: not self-contained, into site/ ---------------------------

dir.create("site", showWarnings = FALSE)
saveWidget(m, "site/map.html", selfcontained = FALSE)
message("Done: site/map.html + site/img/*.png")

if (interactive()) {
  m
}
