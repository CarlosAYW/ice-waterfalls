# =====================================================================
# Icefall thickness model per UID
# - Station: TL/RF (10-min)
# - INCA timeseries: wind + radiation (UU/VV/GL, hourly -> LOCF)
# - Topographic sun + wind-vulnerability LUT
# - Model start always Oct 1 of the current season
# Output: data/ModelRuns/model_uid<UID>.csv
# =====================================================================

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(zoo)
})

TZ_LOCAL <- "Europe/Vienna"

# ----------------------------
# Settings
# ----------------------------
# UID from CLI argument or ENV (GitHub Actions)
args <- commandArgs(trailingOnly = TRUE)
UID_TEST <- if (length(args) >= 1 && nzchar(args[1])) {
  as.integer(args[1])
} else {
  as.integer(Sys.getenv("UID_TEST", "48"))
}
if (!is.finite(UID_TEST)) stop("UID_TEST is not valid.")

MODEL_STEP_MIN <- 10
DT_H   <- MODEL_STEP_MIN / 60
DT_SEC <- MODEL_STEP_MIN * 60
W2MJ_STEP <- DT_SEC / 1e6  # W/m2 -> MJ/m2 per model step

model_end_date_raw <- Sys.getenv("MODEL_END_DATE", "")
END_DATE <- if (nzchar(model_end_date_raw)) as.Date(model_end_date_raw) else Sys.Date()
if (is.na(END_DATE)) stop("MODEL_END_DATE is not a valid date: ", model_end_date_raw)
END_DATE_EXT <- END_DATE

PATH_ASSIGN   <- "data/AWS/icefalls_nearest_station.csv"
PATH_META     <- "data/Koordinaten_Wasserfaelle/eisklettern_links_entries_diff.csv"
PATH_STATIONS <- "data/AWS/stations_all.csv"
DIR_SUN       <- "data/suntime"
PATH_WINDLUT  <- "data/Wind/wind_vulnerability_5deg.csv"
PATH_STRUCTURE <- "data/derived/icefall_structure/icefall_structure_analysis.csv"
PATH_CAP <- "data/CAP/icefall_station_cap.csv"
PATH_INV_DIR <- "data/_cache_inversion"
PATH_INV_RDS <- file.path(PATH_INV_DIR, sprintf("inversion_%s.rds", format(END_DATE, "%Y%m%d")))
PATH_ELEV_DEM <- "data/DEM/DEM_Tirol_INCAgrid_1km_epsg4326.tif"

PATH_INCA_DIR <- "data/inca_nordtirol/point_timeseries"

PATH_OUT      <- sprintf("data/ModelRuns/model_uid%s.csv", UID_TEST)

DEFAULT_SOLAR_TILT_DEG <- 65
MIN_SOLAR_TILT_DEG <- 45
MAX_SOLAR_TILT_DEG <- 85

# ----------------------------
# Helpers
# ----------------------------
season_start_oct <- function(d) {
  d <- as.Date(d)
  y <- year(d); m <- month(d)
  as.Date(sprintf("%d-10-01", ifelse(m >= 10, y, y - 1)))
}
START_DATE <- season_start_oct(END_DATE)

to_num <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

bin5 <- function(dd_deg) {
  ifelse(is.na(dd_deg), NA_real_, (round(dd_deg / 5) * 5) %% 360)
}

parse_time_any <- function(x, tz = TZ_LOCAL) {
  if (inherits(x, "POSIXct")) return(with_tz(x, tz))
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  out <- suppressWarnings(lubridate::ymd_hms(x, tz = tz))
  if (all(is.na(out))) out <- suppressWarnings(lubridate::ymd_hm(x, tz = tz))
  out
}

fill_station_gaps <- function(x, max_gap_steps = as.integer(6 * 60 / MODEL_STEP_MIN)) {
  # Fill only short station outages; never extrapolate long season-start gaps.
  x <- zoo::na.approx(x, na.rm = FALSE, maxgap = max_gap_steps)
  zoo::na.locf(x, na.rm = FALSE, maxgap = max_gap_steps)
}

clamp01 <- function(x) pmin(1, pmax(0, x))

valid_elev_m <- function(x) is.finite(x) & x > -500 & x < 6000

deg2rad <- function(x) x * pi / 180
rad2deg <- function(x) x * 180 / pi

first_finite <- function(...) {
  vals <- c(...)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(NA_real_)
  vals[[1]]
}

haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- 6371
  lat1r <- deg2rad(lat1)
  lat2r <- deg2rad(lat2)
  dlat <- deg2rad(lat2 - lat1)
  dlon <- deg2rad(lon2 - lon1)
  a <- sin(dlat / 2)^2 + cos(lat1r) * cos(lat2r) * sin(dlon / 2)^2
  2 * r * atan2(sqrt(a), sqrt(1 - a))
}

extract_dem_elev_m <- function(lon, lat) {
  if (!file.exists(PATH_ELEV_DEM) || !requireNamespace("terra", quietly = TRUE)) return(NA_real_)
  if (!is.finite(lon) || !is.finite(lat)) return(NA_real_)
  tryCatch({
    dem <- terra::rast(PATH_ELEV_DEM)
    pt <- terra::vect(data.frame(lon = lon, lat = lat), geom = c("lon", "lat"), crs = "EPSG:4326")
    dem_crs <- terra::crs(dem)
    if (nzchar(dem_crs)) pt <- terra::project(pt, dem_crs)
    val <- terra::extract(dem, pt, ID = FALSE)[1, 1]
    val <- to_num(val)
    if (is.finite(val)) val else NA_real_
  }, error = function(e) NA_real_)
}

tz_offset_hours_vec <- function(time_local) {
  z <- format(time_local, "%z")
  sign <- ifelse(substr(z, 1, 1) == "-", -1, 1)
  hh <- suppressWarnings(as.numeric(substr(z, 2, 3)))
  mm <- suppressWarnings(as.numeric(substr(z, 4, 5)))
  off <- sign * (hh + mm / 60)
  off[!is.finite(off)] <- 0
  off
}

read_geosphere_csv <- function(path) {
  x <- tryCatch(readr::read_delim(path, delim = ",", show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (!is.null(x) && ncol(x) > 1) return(x)
  x <- tryCatch(readr::read_delim(path, delim = ";", show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (!is.null(x) && ncol(x) > 1) return(x)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

parse_time_any_utc <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  t <- suppressWarnings(lubridate::ymd_hms(x, tz = "UTC"))
  if (all(is.na(t))) t <- suppressWarnings(lubridate::ymd_hm(x, tz = "UTC"))
  if (all(is.na(t))) t <- suppressWarnings(lubridate::parse_date_time(
    x, orders = c("Ymd HMS", "Ymd HM", "Y-m-d\"T\"H:M:S", "Y-m-d\"T\"H:M"),
    tz = "UTC"
  ))
  t
}

compute_solar_position <- function(time_local, lat_deg, lon_deg) {
  n <- length(time_local)
  if (n == 0 || !is.finite(lat_deg) || !is.finite(lon_deg)) {
    return(list(
      altitude_rad = rep(NA_real_, n),
      altitude_deg = rep(NA_real_, n),
      azimuth_deg = rep(NA_real_, n)
    ))
  }

  doy <- lubridate::yday(time_local)
  hour_local <- lubridate::hour(time_local) +
    lubridate::minute(time_local) / 60 +
    lubridate::second(time_local) / 3600
  gamma <- 2 * pi / 365.2422 * (doy - 1 + (hour_local - 12) / 24)

  eqtime_min <- 229.18 * (
    0.000075 +
      0.001868 * cos(gamma) -
      0.032077 * sin(gamma) -
      0.014615 * cos(2 * gamma) -
      0.040849 * sin(2 * gamma)
  )

  decl <- 0.006918 -
    0.399912 * cos(gamma) +
    0.070257 * sin(gamma) -
    0.006758 * cos(2 * gamma) +
    0.000907 * sin(2 * gamma) -
    0.002697 * cos(3 * gamma) +
    0.00148 * sin(3 * gamma)

  tz_hours <- tz_offset_hours_vec(time_local)
  true_solar_min <- (hour_local * 60 + eqtime_min + 4 * lon_deg - 60 * tz_hours) %% 1440
  hour_angle_deg <- true_solar_min / 4 - 180
  hour_angle_deg[hour_angle_deg < -180] <- hour_angle_deg[hour_angle_deg < -180] + 360

  lat_rad <- deg2rad(lat_deg)
  ha_rad <- deg2rad(hour_angle_deg)

  cos_zen <- sin(lat_rad) * sin(decl) + cos(lat_rad) * cos(decl) * cos(ha_rad)
  cos_zen <- pmin(1, pmax(-1, cos_zen))
  altitude_rad <- asin(cos_zen)

  azimuth_rad <- atan2(
    sin(ha_rad),
    cos(ha_rad) * sin(lat_rad) - tan(decl) * cos(lat_rad)
  )
  azimuth_deg <- (rad2deg(azimuth_rad) + 180) %% 360

  list(
    altitude_rad = altitude_rad,
    altitude_deg = rad2deg(altitude_rad),
    azimuth_deg = azimuth_deg
  )
}

compute_surface_solar_ratio <- function(time_local, lat_deg, lon_deg,
                                        surface_aspect_deg, surface_tilt_deg,
                                        sun_visible,
                                        horiz_floor = 0.15,
                                        max_ratio = 2.6) {
  n <- length(time_local)
  ratio <- ifelse(sun_visible > 0, 1, 0)
  if (n == 0) return(ratio)

  if (!all(is.finite(c(lat_deg, lon_deg, surface_aspect_deg, surface_tilt_deg)))) {
    return(ratio)
  }

  sol <- compute_solar_position(time_local, lat_deg, lon_deg)
  aspect_rad <- deg2rad(surface_aspect_deg %% 360)
  tilt_rad <- deg2rad(pmin(89, pmax(0, surface_tilt_deg)))
  azimuth_rad <- deg2rad(sol$azimuth_deg)

  face_cos <- sin(sol$altitude_rad) * cos(tilt_rad) +
    cos(sol$altitude_rad) * sin(tilt_rad) * cos(azimuth_rad - aspect_rad)
  horiz_cos <- pmax(sin(sol$altitude_rad), horiz_floor)
  raw_ratio <- pmax(0, face_cos) / horiz_cos

  usable <- sun_visible > 0 &
    is.finite(raw_ratio) &
    is.finite(sol$altitude_deg) &
    sol$altitude_deg > 0

  ratio[usable] <- pmin(max_ratio, pmax(1, raw_ratio[usable]))
  ratio
}

# ----------------------------
# INCA timeseries (CSV): UU/VV/GL (hourly)
# ----------------------------
parse_inca_timeseries_csv_file <- function(path, tz_local = TZ_LOCAL) {
  df <- read_geosphere_csv(path)
  if (is.null(df) || nrow(df) == 0) {
    return(tibble(time = as.POSIXct(character(), tz = tz_local), UU = numeric(), VV = numeric(), GL = numeric()))
  }
  
  nms <- names(df)
  nms_low <- tolower(nms)
  time_idx <- which(grepl("time|timestamp|date", nms_low))[1]
  if (is.na(time_idx)) stop("INCA CSV: no time column in ", basename(path))
  time_col <- nms[time_idx]
  
  pick <- function(target) {
    i <- which(nms_low == tolower(target)); if (length(i)) return(nms[i[1]])
    i <- which(grepl(paste0("^", tolower(target), "($|[^a-z0-9])"), nms_low)); if (length(i)) return(nms[i[1]])
    i <- which(grepl(tolower(target), nms_low)); if (length(i)) return(nms[i[1]])
    NA_character_
  }
  
  cUU <- pick("UU"); cVV <- pick("VV"); cGL <- pick("GL")
  
  tibble(
    time_utc = parse_time_any_utc(df[[time_col]]),
    UU = if (!is.na(cUU)) suppressWarnings(as.numeric(df[[cUU]])) else NA_real_,
    VV = if (!is.na(cVV)) suppressWarnings(as.numeric(df[[cVV]])) else NA_real_,
    GL = if (!is.na(cGL)) suppressWarnings(as.numeric(df[[cGL]])) else NA_real_
  ) %>%
    filter(!is.na(time_utc)) %>%
    mutate(time = with_tz(time_utc, tz_local)) %>%
    select(time, UU, VV, GL)
}


download_inca_point_uid_ts <- function(uid, lat, lon, start_date, end_date,
                                       base_dir = PATH_INCA_DIR, verbose = TRUE) {
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  
  cs <- as.Date(start_date); ce <- as.Date(end_date)
  outfile <- file.path(
    base_dir,
    sprintf("inca_uid%s_ts_%s_%s.csv", uid, format(cs, "%Y%m%d"), format(ce, "%Y%m%d"))
  )
  
  if (file.exists(outfile) && file.info(outfile)$size > 0) return(outfile)
  
  base_url <- "https://dataset.api.hub.geosphere.at/v1/timeseries/historical/inca-v1-1h-1km"
  start_time <- sprintf("%sT00:00", format(cs, "%Y-%m-%d"))
  end_time   <- sprintf("%sT23:00", format(ce, "%Y-%m-%d"))
  
  if (verbose) message("INCA TS uid ", uid, ": download ", cs, "..", ce)
  
  resp <- request(base_url) |>
    req_retry(max_tries = 5) |>
    req_error(is_error = function(resp) FALSE) |>
    req_url_query(
      parameters    = "UU,VV,GL",
      start         = start_time,
      end           = end_time,
      lat_lon       = paste0(lat, ",", lon),
      output_format = "csv"
    ) |>
    req_user_agent("icefall-model/1.0 (R httr2)") |>
    req_perform()
  
  if (resp_status(resp) >= 400) {
    msg <- tryCatch(resp_body_string(resp), error = function(e) "")
    stop("INCA download failed (HTTP ", resp_status(resp), "). ", msg)
  }
  
  writeLines(resp_body_string(resp), outfile, useBytes = TRUE)
  Sys.sleep(0.2)
  outfile
}

get_inca_point_hourly <- function(uid, start_date, end_date, lon, lat, path_dir = PATH_INCA_DIR, verbose = TRUE) {
  f <- download_inca_point_uid_ts(uid, lat, lon, start_date, end_date, base_dir = path_dir, verbose = verbose)
  
  parse_inca_timeseries_csv_file(f, tz_local = TZ_LOCAL) %>%
    arrange(time) %>%
    filter(
      time >= as.POSIXct(start_date, tz = TZ_LOCAL),
      time <= as.POSIXct(as.Date(end_date) + 1, tz = TZ_LOCAL) - hours(1)
    ) %>%
    group_by(time) %>%
    summarise(
      UU = mean(UU, na.rm = TRUE),
      VV = mean(VV, na.rm = TRUE),
      GL = mean(GL, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      FF_inca = sqrt(UU^2 + VV^2),
      DD_inca = (atan2(UU, VV) * 180/pi + 180) %% 360,
      GLOW_inca = GL
    ) %>%
    select(time, FF_inca, DD_inca, GLOW_inca)
}

# ----------------------------
# Station weather TL/RF (GeoSphere JSON OR LWD CSV)
# ----------------------------
get_geosphere_station_tlrf <- function(start_date, end_date, station_id) {
  base_url <- "https://dataset.api.hub.geosphere.at/v1/station/historical/klima-v2-10min"
  
  # historical station endpoints expect datetime-like strings (safe choice)
  start_q <- sprintf("%sT00:00", as.character(as.Date(start_date)))
  end_q   <- sprintf("%sT23:50", as.character(as.Date(end_date)))  # 10-min grid
  
  # try both parameter casings (some datasets use lower-case, some upper-case)
  param_tries <- c("tl,rf", "TL,RF")
  
  last_msg <- NULL
  
  for (p in param_tries) {
    resp <- request(base_url) |>
      req_url_query(
        station_ids = as.character(station_id),
        parameters  = p,
        start       = start_q,
        end         = end_q
        # IMPORTANT: do NOT send output_format="json"
      ) |>
      req_user_agent("icefall-model/1.0 (R httr2)") |>
      req_retry(max_tries = 3) |>
      req_error(is_error = function(r) FALSE) |>
      req_perform()
    
    st <- resp_status(resp)
    if (st < 400) {
      dat <- jsonlite::fromJSON(resp_body_string(resp), simplifyVector = FALSE)
      
      ts_raw <- dat[["timestamps"]]
      time_utc <- suppressWarnings(lubridate::ymd_hms(ts_raw, tz = "UTC"))
      if (all(is.na(time_utc))) {
        ts_fix <- ts_raw
        has_colon_tz <- grepl("([+-][0-9]{2}):([0-9]{2})$", ts_fix)
        ts_fix[has_colon_tz] <- paste0(
          substr(ts_fix[has_colon_tz], 1, nchar(ts_fix[has_colon_tz]) - 3),
          substr(ts_fix[has_colon_tz], nchar(ts_fix[has_colon_tz]) - 1, nchar(ts_fix[has_colon_tz]))
        )
        time_utc <- as.POSIXct(strptime(ts_fix, "%Y-%m-%dT%H:%M%z", tz = "UTC"))
      }
      
      feat   <- dat[["features"]][[1]]
      params <- feat[["properties"]][["parameters"]]
      
      pull_anycase <- function(name) {
        # try exact, lower, upper
        cand <- c(name, tolower(name), toupper(name))
        for (nm in cand) {
          p0 <- params[[nm]]
          if (!is.null(p0) && !is.null(p0[["data"]])) {
            lst <- p0[["data"]]
            out <- rep(NA_real_, length(time_utc))
            m <- min(length(lst), length(time_utc))
            out[seq_len(m)] <- vapply(lst[seq_len(m)], function(x) if (is.null(x)) NA_real_ else as.numeric(x), numeric(1))
            return(out)
          }
        }
        rep(NA_real_, length(time_utc))
      }
      
      return(tibble::tibble(
        timestamp = with_tz(time_utc, TZ_LOCAL),
        TL = pull_anycase("tl"),
        RF = pull_anycase("rf")
      ))
    }
    
    # try next casing
    last_msg <- tryCatch(resp_body_string(resp), error = function(e) "")
  }
  
  stop(
    "GeoSphere klima-v2-10min Station-Request failed (tried tl/rf + TL/RF). ",
    "HTTP ", resp_status(resp), ". Response: ", last_msg
  )
}


season_label <- function(date) {
  y <- year(date); m <- month(date)
  ifelse(m >= 10, sprintf("%d_%d", y, y + 1), sprintf("%d_%d", y - 1, y))
}

parse_dt_any <- function(x, tz = TZ_LOCAL) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  out <- suppressWarnings(lubridate::dmy_hms(x, tz = tz))
  if (all(is.na(out))) out <- suppressWarnings(lubridate::dmy_hm(x, tz = tz))
  if (all(is.na(out))) out <- suppressWarnings(lubridate::ymd_hms(x, tz = tz))
  if (all(is.na(out))) out <- suppressWarnings(lubridate::ymd_hm(x, tz = tz))
  out
}

read_lwd_param <- function(station_code, param, season) {
  url <- sprintf("https://wiski.tirol.gv.at/lawine/produkte/ogd/%s/%s_%s_%s.csv",
                 station_code, station_code, param, season)
  
  resp <- tryCatch(request(url) |> req_user_agent("icefall-model/1.0 (R httr2)") |> req_perform(), error = function(e) NULL)
  if (is.null(resp) || resp_status(resp) != 200) return(NULL)
  
  raw <- tryCatch(resp_body_raw(resp), error = function(e) NULL)
  if (is.null(raw) || length(raw) == 0) return(NULL)
  
  txt0 <- tryCatch(rawToChar(raw), error = function(e) NA_character_)
  if (is.na(txt0)) return(NULL)
  txt <- iconv(txt0, from = "", to = "UTF-8", sub = "byte")
  
  tmp <- tryCatch(
    readr::read_delim(file = I(txt), delim = ";", show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (is.null(tmp) || nrow(tmp) == 0) return(NULL)
  
  nms <- names(tmp)
  nms_low <- tolower(iconv(nms, from = "", to = "ASCII//TRANSLIT", sub = ""))
  dt_i <- which(grepl("datetime|date_time|zeitstempel|timestamp|datumzeit|datum_zeit", nms_low))
  dt_col <- if (length(dt_i)) nms[dt_i[1]] else NA_character_
  
  if (!is.na(dt_col)) {
    t <- parse_dt_any(tmp[[dt_col]], tz = TZ_LOCAL)
    val_col <- setdiff(nms, dt_col)[1]
    return(tibble(timestamp = t, value = to_num(tmp[[val_col]])) %>% filter(!is.na(timestamp)))
  }
  
  # fallback: first two cols
  if (ncol(tmp) >= 2) {
    t <- parse_dt_any(tmp[[1]], tz = TZ_LOCAL)
    return(tibble(timestamp = t, value = to_num(tmp[[2]])) %>% filter(!is.na(timestamp)))
  }
  NULL
}

get_lwd_station_tlrf <- function(start_date, end_date, station_code) {
  seasons <- unique(season_label(seq(as.Date(start_date), as.Date(end_date), by = "day")))
  
  tl <- bind_rows(lapply(seasons, function(seas) read_lwd_param(station_code, "LT", seas))) %>%
    mutate(param = "TL")
  rf <- bind_rows(lapply(seasons, function(seas) read_lwd_param(station_code, "LF", seas))) %>%
    mutate(param = "RF")
  
  long <- bind_rows(tl, rf)
  if (nrow(long) == 0) stop("No LWD data for station ", station_code)
  
  wide <- long %>%
    select(timestamp, param, value) %>%
    pivot_wider(names_from = param, values_from = value) %>%
    arrange(timestamp) %>%
    filter(timestamp >= as.POSIXct(start_date, tz = TZ_LOCAL),
           timestamp <  as.POSIXct(end_date + 1, tz = TZ_LOCAL))
  
  for (cc in c("TL","RF")) if (!cc %in% names(wide)) wide[[cc]] <- NA_real_
  wide
}

get_station_tlrf <- function(start_date, end_date, station_id, source) {
  if (identical(source, "GeoSphere")) return(get_geosphere_station_tlrf(start_date, end_date, station_id))
  if (identical(source, "LWD"))      return(get_lwd_station_tlrf(start_date, end_date, station_id))
  stop("Unknown source: ", source)
}

# =====================================================================
# 1) Load inputs (only required columns -> faster)
# =====================================================================
stopifnot(file.exists(PATH_ASSIGN), file.exists(PATH_STATIONS), dir.exists(DIR_SUN), file.exists(PATH_WINDLUT))

assign <- readr::read_csv(PATH_ASSIGN, show_col_types = FALSE, progress = FALSE)
meta_all <- if (file.exists(PATH_META)) {
  readr::read_delim(
    PATH_META, delim = ";",
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE, progress = FALSE,
    col_select = dplyr::any_of(c("uid", "name", "latitude", "longitude", "elevation_dgm5m"))
  )
} else {
  tibble(uid = integer(), name = character(), latitude = numeric(), longitude = numeric(), elevation_dgm5m = numeric())
}
stations_all <- readr::read_csv(
  PATH_STATIONS, show_col_types = FALSE, progress = FALSE,
  col_select = dplyr::any_of(c("station_id", "name", "lon", "lat", "altitude_m", "source", "tl", "rf", "td_lwd"))
)
structure_all <- if (file.exists(PATH_STRUCTURE)) {
  readr::read_csv(
    PATH_STRUCTURE, show_col_types = FALSE, progress = FALSE,
    col_select = dplyr::any_of(c(
      "uid",
      "name",
      "source_aspect_deg",
      "preferred_aspect_deg",
      "slope_mean_deg",
      "slope_p90_deg"
    ))
  )
} else {
  tibble(
    uid = integer(),
    name = character(),
    source_aspect_deg = numeric(),
    preferred_aspect_deg = numeric(),
    slope_mean_deg = numeric(),
    slope_p90_deg = numeric()
  )
}
cap_pairs <- if (file.exists(PATH_CAP)) {
  readr::read_csv(
    PATH_CAP, show_col_types = FALSE, progress = FALSE,
    col_select = dplyr::any_of(c(
      "uid",
      "icefall_cap_potential",
      "station_cap_potential",
      "cap_delta",
      "cap_pair_confidence"
    ))
  )
} else {
  message("CAP cache missing: ", PATH_CAP, " (CAP correction remains 0)")
  tibble(
    uid = integer(),
    icefall_cap_potential = numeric(),
    station_cap_potential = numeric(),
    cap_delta = numeric(),
    cap_pair_confidence = numeric()
  )
}
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

sun_file <- find_sun_file_for_uid(UID_TEST)
if (is.na(sun_file)) {
  message("Warning: missing sun files for UIDs: ", sprintf("%03d", UID_TEST))
  sun_all <- tibble(
    uid = integer(),
    name = character(),
    date = as.Date(character()),
    sunrise_topo = character(),
    sunset_topo = character(),
    sun_hours_topo = numeric()
  )
} else {
  sun_all <- readr::read_csv(
    sun_file, show_col_types = FALSE, progress = FALSE,
    col_select = dplyr::any_of(c("uid","name","date","sunrise_topo","sunset_topo","sun_hours_topo"))
  )
  if (!"sun_hours_topo" %in% names(sun_all)) {
    sun_all$sun_hours_topo <- NA_real_
  }
}
wind_lut <- readr::read_csv(
  PATH_WINDLUT, show_col_types = FALSE, progress = FALSE,
  col_select = dplyr::any_of(c("uid","dir_deg","wind_vuln_0_9"))
)

fallback_station_assignment <- function(uid) {
  ice <- meta_all %>%
    mutate(
      uid_num = suppressWarnings(as.integer(uid)),
      ice_lat = to_num(latitude),
      ice_lon = to_num(longitude),
      icefall_elev_m = to_num(elevation_dgm5m)
    ) %>%
    filter(uid_num == !!uid, is.finite(ice_lat), is.finite(ice_lon)) %>%
    slice(1)

  if (nrow(ice) == 0) return(tibble())

  stations <- stations_all %>%
    mutate(
      station_id = as.character(station_id),
      station_name = as.character(name),
      station_lat = to_num(lat),
      station_lon = to_num(lon),
      station_altitude_m = to_num(altitude_m),
      source = as.character(source)
    ) %>%
    filter(source %in% c("GeoSphere", "LWD"), is.finite(station_lat), is.finite(station_lon))

  if (nrow(stations) == 0) return(tibble())

  stations <- stations %>%
    mutate(dist_km = haversine_km(ice$ice_lat[[1]], ice$ice_lon[[1]], station_lat, station_lon)) %>%
    arrange(dist_km)

  station <- stations %>% slice(1)
  ice_elev_m <- ice$icefall_elev_m[[1]]
  if (!valid_elev_m(ice_elev_m)) {
    ice_elev_m <- extract_dem_elev_m(ice$ice_lon[[1]], ice$ice_lat[[1]])
  }
  if (!valid_elev_m(ice_elev_m)) {
    ice_elev_m <- station$station_altitude_m[[1]]
  }
  elev_diff_m <- ice_elev_m - station$station_altitude_m[[1]]

  tibble(
    uid = uid,
    icefall_name = as.character(ice$name[[1]]),
    ice_lat = ice$ice_lat[[1]],
    ice_lon = ice$ice_lon[[1]],
    station_id = station$station_id[[1]],
    station_name = station$station_name[[1]],
    source = station$source[[1]],
    station_lat = station$station_lat[[1]],
    station_lon = station$station_lon[[1]],
    dist_km = station$dist_km[[1]],
    icefall_elev_m = ice_elev_m,
    elev_diff_m = elev_diff_m,
    topo_pos_diff = NA_real_,
    path_barrier_m = NA_real_
  )
}

row_uid <- assign %>% filter(uid == UID_TEST) %>% slice(1)
if (nrow(row_uid) == 0) {
  message("UID ", UID_TEST, " not found in station assignment table; using nearest-station fallback.")
  row_uid <- fallback_station_assignment(UID_TEST)
}
if (nrow(row_uid) == 0) stop("uid not found: ", UID_TEST)

station_id <- as.character(row_uid$station_id)
source     <- as.character(row_uid$source)
dist_km    <- to_num(row_uid$dist_km)
dz_m       <- to_num(row_uid$elev_diff_m)  # ice - station
ice_lon    <- to_num(row_uid$ice_lon)
ice_lat    <- to_num(row_uid$ice_lat)
path_barrier_m_uid <- if ("path_barrier_m" %in% names(row_uid)) to_num(row_uid$path_barrier_m) else NA_real_
topo_pos_diff_uid <- if ("topo_pos_diff" %in% names(row_uid)) to_num(row_uid$topo_pos_diff) else NA_real_
if (!is.finite(path_barrier_m_uid)) path_barrier_m_uid <- 80
if (!is.finite(topo_pos_diff_uid)) topo_pos_diff_uid <- 0.16

cap_uid <- cap_pairs %>% filter(uid == UID_TEST) %>% slice(1)
cap_num <- function(df, col, default = NA_real_) {
  if (nrow(df) == 1 && col %in% names(df)) return(to_num(df[[col]]))
  default
}
cap_icefall_uid <- cap_num(cap_uid, "icefall_cap_potential")
cap_station_uid <- cap_num(cap_uid, "station_cap_potential")
cap_delta_uid <- cap_num(cap_uid, "cap_delta", 0)
cap_pair_confidence_uid <- cap_num(cap_uid, "cap_pair_confidence", 0)
cap_delta_uid <- ifelse(is.finite(cap_delta_uid), pmax(-1, pmin(1, cap_delta_uid)), 0)
cap_pair_confidence_uid <- ifelse(is.finite(cap_pair_confidence_uid), clamp01(cap_pair_confidence_uid), 0)

ice_name <- NA_character_

if ("icefall_name" %in% names(row_uid)) ice_name <- row_uid$icefall_name
if (is.null(ice_name) || is.na(ice_name) || !nzchar(ice_name)) {
  if ("ice_name" %in% names(row_uid)) ice_name <- row_uid$ice_name
}
if (is.null(ice_name) || is.na(ice_name) || !nzchar(ice_name)) ice_name <- row_uid$name
if (is.null(ice_name) || is.na(ice_name) || !nzchar(ice_name)) ice_name <- "Icefall"

ice_alt_m <- if ("icefall_elev_m" %in% names(row_uid)) to_num(row_uid$icefall_elev_m) else NA_real_
ice_fallheight_m <- if ("icefall_height_m" %in% names(row_uid)) to_num(row_uid$icefall_height_m) else NA_real_

st_meta <- stations_all %>% filter(as.character(station_id) == .env$station_id) %>% slice(1)
z_aws <- if (nrow(st_meta) == 1) to_num(st_meta$altitude_m) else NA_real_

structure_uid <- structure_all %>% filter(uid == UID_TEST) %>% slice(1)
ice_aspect_deg_uid <- if (nrow(structure_uid) == 1) {
  first_finite(
    to_num(structure_uid$preferred_aspect_deg),
    to_num(structure_uid$source_aspect_deg)
  )
} else {
  NA_real_
}
ice_slope_mean_deg_uid <- if (nrow(structure_uid) == 1) to_num(structure_uid$slope_mean_deg) else NA_real_
ice_slope_p90_deg_uid <- if (nrow(structure_uid) == 1) to_num(structure_uid$slope_p90_deg) else NA_real_
ice_tilt_seed_deg_uid <- first_finite(
  if (is.finite(ice_slope_mean_deg_uid) && is.finite(ice_slope_p90_deg_uid)) {
    0.35 * ice_slope_mean_deg_uid + 0.65 * ice_slope_p90_deg_uid
  } else {
    NA_real_
  },
  ice_slope_p90_deg_uid,
  ice_slope_mean_deg_uid,
  DEFAULT_SOLAR_TILT_DEG
)
ice_tilt_deg_uid <- pmin(MAX_SOLAR_TILT_DEG, pmax(MIN_SOLAR_TILT_DEG, ice_tilt_seed_deg_uid))

# =====================================================================
# 2) Sun + Wind LUT
# =====================================================================
sun_uid <- sun_all %>%
  filter(uid == UID_TEST) %>%
  mutate(
    date = as.Date(date),
    sunrise_topo = parse_time_any(sunrise_topo, tz = TZ_LOCAL),
    sunset_topo  = parse_time_any(sunset_topo,  tz = TZ_LOCAL),
    sun_hours_topo = to_num(sun_hours_topo)
  ) %>%
  select(date, sunrise_topo, sunset_topo, sun_hours_topo) %>%
  distinct(date, .keep_all = TRUE) %>%
  filter(date >= START_DATE, date <= END_DATE_EXT)

if (all(is.na(sun_uid$sun_hours_topo))) {
  sun_uid <- sun_uid %>%
    mutate(
      sun_hours_topo = as.numeric(difftime(sunset_topo, sunrise_topo, units = "hours"))
    )
}

if (nrow(sun_uid) == 0) message("Warning: no sun data for uid ", UID_TEST, " (topo_sun_fac remains 0)")

wind_uid <- wind_lut %>%
  filter(uid == UID_TEST) %>%
  mutate(dir_deg = as.numeric(dir_deg), wind_vuln_0_9 = as.integer(wind_vuln_0_9)) %>%
  select(dir_deg, wind_vuln_0_9)

if (nrow(wind_uid) == 0) {
  message("No wind LUT for uid ", UID_TEST, " (wind_vuln_0_9=9 for all directions)")
  wind_uid <- tibble(
    dir_deg = seq(0, 355, by = 5),
    wind_vuln_0_9 = 9L
  )
}

# =====================================================================
# 3) Station TL/RF (10-min) -> 10-min Raster + Fill
# =====================================================================
wx10 <- get_station_tlrf(START_DATE, END_DATE, station_id, source) %>%
  mutate(
    timestamp = with_tz(as.POSIXct(timestamp, tz = TZ_LOCAL), TZ_LOCAL),
    TL = to_num(TL),
    RF = to_num(RF)
  )

if (all(is.na(wx10$TL))) stop("TL is completely missing.")
if (all(is.na(wx10$RF))) stop("RF is completely missing.")

step_str <- paste0(MODEL_STEP_MIN, " mins")
model_time <- seq(
  as.POSIXct(START_DATE, tz = TZ_LOCAL),
  as.POSIXct(END_DATE + 1, tz = TZ_LOCAL) - lubridate::minutes(MODEL_STEP_MIN),
  by = step_str
)

wx <- wx10 %>%
  mutate(time = floor_date(timestamp, unit = step_str)) %>%
  group_by(time) %>%
  summarise(
    TL = if (all(is.na(TL))) NA_real_ else mean(TL, na.rm = TRUE),
    RF = if (all(is.na(RF))) NA_real_ else mean(RF, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::complete(time = model_time) %>%
  arrange(time) %>%
  mutate(
    TL = fill_station_gaps(TL),
    RF = fill_station_gaps(RF)
  )

# =====================================================================
# 4) INCA override: FF/DD/GLOW (hourly -> 10-min LOCF)
# =====================================================================
inca_hr <- get_inca_point_hourly(
  uid = UID_TEST,
  start_date = START_DATE,
  end_date   = END_DATE,
  lon = ice_lon, lat = ice_lat,
  path_dir = PATH_INCA_DIR,
  verbose = TRUE
)

wx <- wx %>%
  left_join(inca_hr, by = "time") %>%
  arrange(time) %>%
  mutate(
    FF_inca   = zoo::na.locf(FF_inca,   na.rm = FALSE),
    DD_inca   = zoo::na.locf(DD_inca,   na.rm = FALSE),
    GLOW_inca = zoo::na.locf(GLOW_inca, na.rm = FALSE),
    FF_inca   = zoo::na.locf(FF_inca,   na.rm = FALSE, fromLast = TRUE),
    DD_inca   = zoo::na.locf(DD_inca,   na.rm = FALSE, fromLast = TRUE),
    GLOW_inca = zoo::na.locf(GLOW_inca, na.rm = FALSE, fromLast = TRUE),
    FF   = coalesce(FF_inca, 0),
    DD   = coalesce(DD_inca, 0),
    GLOW = coalesce(GLOW_inca, 0)
  ) %>%
  select(-FF_inca, -DD_inca, -GLOW_inca) %>%
  filter(
    time >= as.POSIXct(START_DATE, tz = TZ_LOCAL),
    time <  as.POSIXct(END_DATE + 1, tz = TZ_LOCAL)
  )

# --- Inversion join (global, same time base for all UIDs) ---
if (file.exists(PATH_INV_RDS)) {
  inv <- readRDS(PATH_INV_RDS) %>%
    mutate(time = as.POSIXct(time, tz = TZ_LOCAL))
  wx <- wx %>%
    left_join(inv %>% select(time, inv_active, inv_score_C, inv_class,
                             inv_grad_max_C_per_100m,
                             grad01_K_per_m, grad12_K_per_m, grad02_K_per_m),
              by = "time") %>%
    mutate(
      inv_active = ifelse(is.na(inv_active), FALSE, inv_active),
      inv_score_C  = ifelse(is.na(inv_score_C),  0,     inv_score_C),
      inv_class  = ifelse(is.na(inv_class),  "none", inv_class)
    )
} else {
  wx <- wx %>% mutate(
    inv_active = FALSE,
    inv_score_C = 0,
    inv_class = "none",
    inv_grad_max_C_per_100m = NA_real_,
    grad01_K_per_m = NA_real_,
    grad12_K_per_m = NA_real_,
    grad02_K_per_m = NA_real_
  )
  message("Warning: inversion cache missing: ", PATH_INV_RDS, " (inv_active=FALSE)")
}

# =====================================================================
# 6) Join topo sun + wind vulnerability
# =====================================================================
wx <- wx %>%
  mutate(date = as.Date(time)) %>%
  left_join(sun_uid, by = "date") %>%
  mutate(
    topo_sun_fac = ifelse(
      !is.na(sunrise_topo) & !is.na(sunset_topo) & time >= sunrise_topo & time < sunset_topo,
      1, 0
    ),
    dd_bin = bin5(DD)
  ) %>%
  left_join(wind_uid, by = c("dd_bin" = "dir_deg")) %>%
  mutate(
    wind_vuln_0_9 = ifelse(is.na(wind_vuln_0_9), 9L, wind_vuln_0_9),
    wind_vuln = pmin(1, pmax(0, wind_vuln_0_9 / 9))
  )

wx$solar_incidence_ratio <- compute_surface_solar_ratio(
  time_local = wx$time,
  lat_deg = ice_lat,
  lon_deg = ice_lon,
  surface_aspect_deg = ice_aspect_deg_uid,
  surface_tilt_deg = ice_tilt_deg_uid,
  sun_visible = wx$topo_sun_fac
)

# =====================================================================
# 7) Ice model
# =====================================================================
ice_params <- list(
  albedo = 0.50,              # Shortwave albedo of the ice surface [-];
  # fraction of incoming solar radiation reflected by the ice.
  # Higher values reduce absorbed radiation and therefore reduce melt.

  solar_boost_strength = 1.35,# Extra shortwave boost for well-exposed wall geometry [-];
  # south/east/west walls absorb more radiation than a horizontal surface
  # when they are directly illuminated.

  solar_core_boost_strength = 0.55,
  # converts the geometric solar boost into stronger core/reserve melt.
  # kept lower than the surface-radiation boost to avoid overreacting.

  solar_load_cap = 3.0,       # Upper cap for solar melt amplification [-];
  # prevents unrealistically large boosts at very low sun angles.

  solar_core_cap = 1.8,       # Upper cap for core/reseve solar forcing [-];
  # keeps sun-driven deep melt bounded even on very exposed lines.
  
  Hmax_m = 0.90,              # Asymptotic maximum ice thickness [m];
  # used as the upper thickness scale for growth saturation.
  
  H0_m = 0.20,                # Exponential insulation / damping scale [m];
  # controls how quickly further ice growth is reduced as total ice
  # thickness increases. Smaller values = stronger early damping.
  
  Hsat_lambda_m = 0.15,       # Exponential saturation softness parameter [m];
  # controls how smoothly growth approaches zero near Hmax_m.
  # Smaller values = sharper cutoff near Hmax_m,
  # larger values = smoother, more gradual saturation.
  
  core_seed_mm = 4,           # Baseline initial core-ice thickness at model start [mm];
  # represents a small amount of pre-existing dense/supportive inner ice.
  
  core_seed_boost_mm = 55,    # Additional core seed for strong cold-retention sites [mm];
  # large vertical station offsets justify more persistent pre-existing core ice.
  
  surface_seed_mm = 16,       # Initial surface-ice thickness at model start [mm];
  # represents pre-existing outer / recently formed surface ice.
  
  core_share_min = 0.02,      # Baseline minimum fraction of new growth routed to core ice [-];
  # low-retention sites keep less of their new ice in the persistent core layer.
  
  core_share_min_boost = 0.22,# Extra minimum core-growth share for strong cold-retention sites [-];
  # cold-retention sites convert more new ice into long-lived core ice.
  
  core_share_max = 0.10,      # Baseline maximum fraction of new growth routed to core ice [-];
  # limits core routing for low-retention sites.
  
  core_share_max_boost = 0.48,# Extra maximum core-growth share for strong cold-retention sites [-];
  # allows deeply retained sites to build a much more persistent core.
  
  melt_damping_max = 0.05,    # Baseline reduction of melt for thick ice [-];
  # low-retention sites receive only moderate melt protection.
  
  melt_damping_boost = 0.50,  # Additional melt protection for strong cold-retention sites [-];
  # thick ice in strongly retained sites stays shielded much longer.
  
  melt_damping_scale_m = 0.45,# Thickness scale [m] over which melt damping increases;
  # around this thickness, the protection effect of thicker ice
  # becomes substantial.
  
  core_melt_base = 0.48,      # Baseline fraction of residual melt that can affect the core [-];
  # low-retention sites expose the core more readily once surface ice is gone.
  
  core_melt_base_drop = 0.32, # Reduction of base core-melt sensitivity for retained sites [-];
  # strongly retained sites keep the core better insulated from residual melt.
  
  core_melt_warm = 0.30,      # Baseline extra core-melt sensitivity to positive air temperature [-];
  # warm air erodes the core more aggressively at low-retention sites.
  
  core_melt_warm_drop = 0.12, # Reduction of warm-air core-melt sensitivity for retained sites [-];
  # retained sites respond less strongly to the same positive temperatures.
  
  core_melt_sun = 0.24,       # Baseline extra core-melt sensitivity when the icefall is sunlit [-];
  # direct topographic sun exposure increases core vulnerability.
  
  core_melt_sun_drop = 0.10   # Reduction of sun-driven core melt for retained sites [-];
  # retained sites lose less core ice under the same sun exposure.
)

coef <- list(
  lapse_K_per_m      = 0.0065,# Environmental lapse rate [K m^-1];
  # used to vertically adjust air temperature from station/grid
  # elevation to icefall elevation.
  
  growth_mm_per_C_h  = 0.50,  # Ice growth coefficient [mm per °C per hour];
  # converts freezing degree hours into potential ice growth.
  
  melt_mm_per_C_h    = 0.60,  # Melt coefficient [mm per °C per hour];
  # converts positive degree hours into potential melt.
  
  rad_melt_mm_per_MJ = 0.35,  # Radiation melt coefficient [mm per MJ m^-2];
  # converts absorbed shortwave energy into additional melt.
  
  k_wind             = 0.06,  # Wind enhancement factor [- per m s^-1];
  # scales growth/melt intensity with wind speed
  
  k_dry              = 0.25,  # Dryness enhancement factor [-];
  # increases growth potential under dry-air conditions
  # through a relative-humidity-based multiplier.
  
  wind_cap_ms        = 15,     # Maximum wind speed considered by the model [m s^-1];
  # caps the wind effect to avoid unrealistically large responses
  # at very high wind speeds.

  cap_max_adjust_C   = 6.5,    # Maximum relative CAP temperature correction [K];
  # applied as station-relative cold-air-pooling potential.
  # Anchored in Rauchoecker et al. (2024, QJRMS 150, 1243-1266), who measured
  # 5-9 K between valley-floor and slope stations in a 200-300 m deep Alpine
  # valley (SEECAP, Seefeld/Tyrol); 6.5 K is the midpoint of that range.

  cap_wind_full_ms   = 1.0,    # Wind speed up to which the CAP is undisturbed [m s^-1];
  # Rauchoecker et al. (2024) observed warming inside the pool already above
  # ~1 m s^-1.

  cap_wind_shutdown_ms = 2.0,  # Wind speed where CAP cooling is mixed out [m s^-1];
  # Rauchoecker et al. (2024): the pool was eroded completely above ~2 m s^-1
  # (foehn breakthrough).

  cap_radiation_shutdown_Wm2 = 180, # Direct-radiation threshold reducing CAP [-];
  # direct sun weakens near-surface pooling during the day.

  cap_tau_build_h    = 2.9,    # Cooling time constant of the pool build-up [h];
  # Rauchoecker et al. (2024) report 2.9 h for Seefeld (Gruenloch 3.2 h,
  # Peter Sinks 3.1 h) with an initial cooling rate of 3.0 K h^-1. The pool
  # therefore needs several hours to reach full strength after onset.

  cap_no_inversion_factor = 0.35 # Residual CAP strength when no inversion is detected [-].
  # keeps weak nocturnal/local pooling possible even without the global inversion flag.
)

# Hours since the regional inversion switched on, used for the exponential
# build-up of the cold-air pool (cap_tau_build_h). Counts consecutive
# inv_active steps and resets whenever the inversion is off.
.cap_on <- ifelse(is.na(wx$inv_active), FALSE, as.logical(wx$inv_active))
wx$cap_since_onset_h <- stats::ave(
  as.integer(.cap_on), cumsum(!.cap_on), FUN = cumsum
) * MODEL_STEP_MIN / 60

wx <- wx %>%
  mutate(
    dz_raw = dz_m,
    
    # Apply the physical profile only when inversion is active and station elevation is known.
    use_prof = inv_active & is.finite(z_aws) & is.finite(dz_m) &
      is.finite(grad01_K_per_m) & is.finite(grad12_K_per_m),
    
    # Target elevation from station elevation + dz (also works if ice_alt_m is missing).
    z_target_m = z_aws + dz_m,
    
    # Piecewise dT over two layers (below/above Z1=1935 m).
    dT_prof = if_else(
      use_prof,
      {
        z1_m <- 1935
        lo <- pmin(z_aws, z_target_m)
        hi <- pmax(z_aws, z_target_m)
        
        len_low  <- pmax(0, pmin(hi, z1_m) - lo)
        len_high <- pmax(0, hi - pmax(lo, z1_m))
        
        sgn <- if_else(z_target_m >= z_aws, 1, -1)
        sgn * (grad01_K_per_m * len_low + grad12_K_per_m * len_high)
      },
      NA_real_
    ),
    
    # Temperature at icefall elevation:
    # - Default: constant lapse
    # - If inv_active: physical profile (also for dz > 0 and dz < 0)
    TLz_raw = TL - coef$lapse_K_per_m * dz_raw,
    TLz_base = if_else(use_prof, TL + dT_prof, TLz_raw),
    FF_eff = pmin(coef$wind_cap_ms, pmax(0, FF)),
    GLOW = if_else(is.finite(GLOW), GLOW, 0),
    weather_ok = is.finite(TL),

    # Station-relative cold-air pooling correction:
    # positive cap_delta_uid cools the route relative to the station;
    # negative values warm it when the station is the stronger cold-air pool.
    # Wind gating: full strength up to cap_wind_full_ms, linearly to zero at
    # cap_wind_shutdown_ms (Rauchoecker et al. 2024).
    cap_wind_fac = clamp01(
      (coef$cap_wind_shutdown_ms - FF_eff) /
        (coef$cap_wind_shutdown_ms - coef$cap_wind_full_ms)
    ),
    cap_sun_fac = if_else(
      topo_sun_fac > 0,
      clamp01(1 - GLOW / coef$cap_radiation_shutdown_Wm2),
      1
    ),
    cap_inv_fac = if_else(inv_active, 1, coef$cap_no_inversion_factor),
    # Exponential build-up after inversion onset (tau = 2.9 h, Seefeld).
    # Only active while the regional inversion is on; otherwise the residual
    # cap_no_inversion_factor already represents the weak local pooling.
    cap_ramp_fac = if_else(
      inv_active,
      clamp01(1 - exp(-cap_since_onset_h / coef$cap_tau_build_h)),
      1
    ),
    cap_stability_fac = cap_wind_fac * cap_sun_fac * cap_inv_fac * cap_ramp_fac,
    cap_temp_adjust_C = coef$cap_max_adjust_C * cap_delta_uid *
      cap_pair_confidence_uid * cap_stability_fac,
    TLz = TLz_base - cap_temp_adjust_C,

    FDH = if_else(is.finite(TLz), pmax(0, -TLz), 0),
    PDH = if_else(is.finite(TLz), pmax(0,  TLz), 0),

    solar_load_fac = if_else(
      topo_sun_fac > 0,
      pmin(
        ice_params$solar_load_cap,
        1 + ice_params$solar_boost_strength * pmax(0, solar_incidence_ratio - 1)
      ),
      0
    ),
    solar_core_fac = if_else(
      topo_sun_fac > 0,
      pmin(
        ice_params$solar_core_cap,
        1 + ice_params$solar_core_boost_strength * pmax(0, solar_incidence_ratio - 1)
      ),
      0
    ),
    SW_MJ_step = GLOW * W2MJ_STEP * solar_load_fac * (1 - ice_params$albedo),
    
    wind_fac = 1 + coef$k_wind * FF_eff * wind_vuln,
    dry_fac  = if_else(is.finite(RF), 1 + coef$k_dry * pmax(0, 1 - RF / 100), 1),
    
    base_growth_mm_step = if_else(
      weather_ok,
      coef$growth_mm_per_C_h * FDH * DT_H * wind_fac * dry_fac,
      0
    ),
    base_melt_mm_step = if_else(
      weather_ok,
      coef$melt_mm_per_C_h * PDH * DT_H * wind_fac +
        coef$rad_melt_mm_per_MJ * SW_MJ_step,
      0
    )
  ) %>%
  mutate(
    TLz_72h_step = zoo::rollapplyr(
      TLz, width = as.integer(72 * 60 / MODEL_STEP_MIN),
      FUN = function(x) {
        x <- x[is.finite(x)]
        if (length(x) == 0) NA_real_ else mean(x)
      },
      fill = NA_real_, partial = TRUE
    )
  )

surface_mm <- numeric(nrow(wx))
core_mm <- numeric(nrow(wx))
retention_dz_fac_uid <- clamp01((abs(dz_m) - 150) / 300)
retention_barrier_fac_uid <- clamp01((path_barrier_m_uid - 80) / 160)
retention_topo_fac_uid <- clamp01((0.16 - topo_pos_diff_uid) / 0.10)
retention_fac_uid <- retention_dz_fac_uid * retention_barrier_fac_uid * retention_topo_fac_uid
exposure_barrier_fac_uid <- 1 - clamp01((path_barrier_m_uid - 40) / 120)
exposure_topo_fac_uid <- clamp01((topo_pos_diff_uid - 0.18) / 0.12)
exposure_fac_uid <- exposure_barrier_fac_uid * exposure_topo_fac_uid
core_reserve_mm <- numeric(nrow(wx))

surface_mm[1] <- ice_params$surface_seed_mm
core_mm[1] <- ice_params$core_seed_mm + ice_params$core_seed_boost_mm * retention_fac_uid

for (i in 2:nrow(wx)) {
  Hprev_m <- (surface_mm[i-1] + core_mm[i-1]) / 1000
  iso <- exp(-Hprev_m / ice_params$H0_m)
  
  gap_m <- pmax(0, ice_params$Hmax_m - Hprev_m)
  cap_exp <- (1 - exp(-gap_m / ice_params$Hsat_lambda_m)) /
    (1 - exp(-ice_params$Hmax_m / ice_params$Hsat_lambda_m))
  
  growth_total <- wx$base_growth_mm_step[i] * iso * cap_exp

  core_share_min_eff <- ice_params$core_share_min +
    ice_params$core_share_min_boost * retention_fac_uid
  core_share_max_eff <- ice_params$core_share_max +
    ice_params$core_share_max_boost * retention_fac_uid
  core_share <- core_share_min_eff +
    0.14 * pmin(1, wx$FDH[i] / 4) +
    0.12 * pmin(1, Hprev_m / 0.30)
  core_share <- pmin(core_share_max_eff, pmax(core_share_min_eff, core_share))

  growth_core <- growth_total * core_share
  growth_surface <- growth_total - growth_core

  surface_pre_melt <- surface_mm[i-1] + growth_surface
  core_pre_melt <- core_mm[i-1] + growth_core

  melt_damping_eff <- ice_params$melt_damping_max +
    ice_params$melt_damping_boost * retention_fac_uid
  melt_scale <- 1 - melt_damping_eff * pmin(1, Hprev_m / ice_params$melt_damping_scale_m)
  melt_total <- wx$base_melt_mm_step[i] * melt_scale
  spring_exposure_fac <- clamp01((wx$TLz_72h_step[i] + 5) / 5)
  if (!is.finite(spring_exposure_fac)) spring_exposure_fac <- 0
  melt_total <- melt_total * (1 + 0.8 * exposure_fac_uid * spring_exposure_fac)

  melt_surface <- min(surface_pre_melt, melt_total)
  melt_left <- pmax(0, melt_total - melt_surface)

  core_melt_fac <- (ice_params$core_melt_base - ice_params$core_melt_base_drop * retention_fac_uid) +
    (ice_params$core_melt_warm - ice_params$core_melt_warm_drop * retention_fac_uid) * pmin(1, wx$PDH[i] / 4) +
    (ice_params$core_melt_sun - ice_params$core_melt_sun_drop * retention_fac_uid) * wx$solar_core_fac[i]
  core_melt_fac <- pmin(0.85, pmax(0.20, core_melt_fac))
  melt_core <- min(core_pre_melt, melt_left * core_melt_fac)

  reserve_gain_mm <- 0.30 * growth_core * retention_fac_uid
  reserve_temp_loss <- pmax(0, wx$TLz_72h_step[i] + 2)
  if (!is.finite(reserve_temp_loss)) reserve_temp_loss <- 0
  reserve_loss_mm <- (0.08 * wx$PDH[i] * DT_H + 0.05 * wx$solar_core_fac[i] + 0.02 * reserve_temp_loss) *
    retention_fac_uid
  reserve_cap_mm <- 0.28 * core_pre_melt
  core_reserve_mm[i] <- min(
    reserve_cap_mm,
    max(0, core_reserve_mm[i-1] + reserve_gain_mm - reserve_loss_mm)
  )

  surface_mm[i] <- max(0, surface_pre_melt - melt_surface)
  core_mm[i] <- max(core_reserve_mm[i], core_pre_melt - melt_core)
}

mod <- wx %>%
  mutate(
    thickness_m = (surface_mm + core_mm) / 1000,
    surface_ice_m = surface_mm / 1000,
    core_ice_m = core_mm / 1000,
    station_id = station_id,
    source = source,
    dist_km = dist_km,
    dz_m = dz_m,
    cap_icefall = cap_icefall_uid,
    cap_station = cap_station_uid,
    cap_delta = cap_delta_uid,
    cap_pair_confidence = cap_pair_confidence_uid
  )

# =====================================================================
# 8) Save
# =====================================================================
mod <- mod %>% mutate(date = as.Date(time))

dir.create(dirname(PATH_OUT), showWarnings = FALSE, recursive = TRUE)
write_csv(mod, PATH_OUT)

# =====================================================================
# 9) Plot historical ice thickness
# =====================================================================
x_min <- as.POSIXct(START_DATE, tz = TZ_LOCAL)
x_max <- max(mod$time, na.rm = TRUE)

plt <- ggplot(mod, aes(time, thickness_m)) +
  geom_line(color = "black", linewidth = 0.9, na.rm = TRUE) +
  coord_cartesian(xlim = c(x_min, x_max)) +
  scale_x_datetime(
    date_breaks = "1 month",
    date_labels = "%b",
    timezone = TZ_LOCAL,
    guide = guide_axis(check.overlap = TRUE)
  ) +
  scale_y_continuous(name = "Ice thickness (m)") +
  labs(
    title = paste0("Modeled ice thickness - ", ice_name, " (UID ", sprintf("%03d", UID_TEST), ")"),
    subtitle = paste(
      c(
        if (!is.na(ice_fallheight_m)) paste0("Icefall height: ", round(ice_fallheight_m, 0), " m"),
        if (!is.na(ice_alt_m)) paste0("Elevation: ", round(ice_alt_m, 0), " m"),
        paste0("Station: ", station_id, " (", source, ")"),
        paste0("dist ", round(dist_km, 2), " km"),
        paste0("dz ", round(dz_m, 0), " m"),
        if (cap_pair_confidence_uid > 0) paste0("CAP delta ", round(cap_delta_uid, 2))
      ),
      collapse = " | "
    ),
    caption = paste0("10-min model (dt=", MODEL_STEP_MIN, " min): FDH/PDH + CAP + SW(toposun) + Wind(vuln) + Dryness + Saturation"),
    x = "Time"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 9, lineheight = 0.95),
    plot.title = element_text(face = "bold")
  )

# =====================================================================
# Export for the web: PNG to site/plots/
# =====================================================================
dir.create("site/plots", recursive = TRUE, showWarnings = FALSE)

plot_file <- sprintf("site/plots/uid_%03d.png", UID_TEST)


ggsave(
  filename = plot_file,
  plot     = plt,
  width    = 14,
  height   = 5,
  units    = "in",
  dpi      = 200,
  bg       = "white"
)

message("Plot written: ", plot_file)
