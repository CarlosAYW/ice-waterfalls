# ============================================================
# Build ONE global inversion time series (not per-UID) and cache it.
# - GeoSphere station: 38 (Imst) 10-min TL
# - LWD stations: IMUT2 (Muttekopfhütte) TL, IMUT1 (Vorderer Maldonkopf) TL
# Output:
#   data/_cache_inversion/inversion_YYYYMMDD.rds
# ============================================================

suppressPackageStartupMessages({
  library(httr2)
  library(httr)
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(readr)
  library(zoo)
  library(stringr)
})

TZ_LOCAL <- "Europe/Vienna"
MODEL_STEP_MIN <- 10
step_str <- paste0(MODEL_STEP_MIN, " mins")

# ----------------------------
# Time range (as in the model)
# ----------------------------
NOW_LOCAL <- with_tz(Sys.time(), TZ_LOCAL)
model_end_date_raw <- Sys.getenv("MODEL_END_DATE", "")
END_DATE <- if (nzchar(model_end_date_raw)) as.Date(model_end_date_raw) else as.Date(NOW_LOCAL)
if (is.na(END_DATE)) stop("MODEL_END_DATE is not a valid date: ", model_end_date_raw)

season_start_oct <- function(d) {
  d <- as.Date(d)
  y <- year(d); m <- month(d)
  as.Date(sprintf("%d-10-01", ifelse(m >= 10, y, y - 1)))
}
# Full season (Oct 1 through today) so inv_active exists for every 10-min timestamp
START_DATE <- season_start_oct(END_DATE)

# ----------------------------
# Helpers
# ----------------------------

Z0_m <-  860  # Imst
Z1_m <- 1935  # Muttekopfhütte
Z2_m <- 2580  # Vorderer Maldonkopf


`%||%` <- function(a, b) if (!is.null(a)) a else b

to_num <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

# --- robust datetime parse without lubridate::parse_date_time (avoids regex issues) ---
parse_dt_safe <- function(x, tz = TZ_LOCAL) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  
  # try common formats seen in LWD/GeoSphere CSVs
  fmts <- c(
    "%Y-%m-%dT%H:%M:%OSZ",
    "%Y-%m-%dT%H:%M:%OS",
    "%Y-%m-%d %H:%M:%OS",
    "%d.%m.%Y %H:%M:%OS",
    "%d.%m.%Y %H:%M",
    "%Y-%m-%d %H:%M",
    "%Y-%m-%d"
  )
  
  out <- rep(as.POSIXct(NA), length(x))
  for (f in fmts) {
    miss <- is.na(out) & !is.na(x)
    if (!any(miss)) break
    cand <- suppressWarnings(as.POSIXct(x[miss], tz = tz, format = f))
    out[miss] <- cand
  }
  
  # if still NA: try lubridate ymd_hms/dmy_hms (these usually work without heavy regex)
  miss <- is.na(out) & !is.na(x)
  if (any(miss)) {
    cand <- suppressWarnings(lubridate::ymd_hms(x[miss], tz = tz))
    out[miss] <- cand
  }
  miss <- is.na(out) & !is.na(x)
  if (any(miss)) {
    cand <- suppressWarnings(lubridate::dmy_hms(x[miss], tz = tz))
    out[miss] <- cand
  }
  
  out
}

# ----------------------------
# GeoSphere TL (CSV output to avoid JSON timestamp parsing)
# ----------------------------
get_geosphere_tl <- function(start_date, end_date, station_id) {
  # dataset: klima-v2-10min (per GeoSphere changelog)
  base_url <- "https://dataset.api.hub.geosphere.at/v1/station/historical/klima-v2-10min"
  
  # Use full-day window; we can cut later.
  start_q <- sprintf("%s", as.character(as.Date(start_date)))
  end_q   <- sprintf("%s", as.character(as.Date(end_date + 1)))  # end is exclusive (next day)
  
  req <- request(base_url) |>
    req_url_query(
      station_ids    = as.character(station_id),
      parameters     = "tl",
      start          = start_q,
      end            = end_q,
      output_format  = "csv"
    ) |>
    req_user_agent("icefall-model/1.0 (R httr2)") |>
    req_retry(max_tries = 3)
  
  resp <- req_perform(req)
  if (resp_status(resp) >= 400) {
    stop("GeoSphere TL failed HTTP ", resp_status(resp), ": ", resp_body_string(resp))
  }
  
  txt <- resp_body_string(resp)
  df <- readr::read_csv(I(txt), show_col_types = FALSE, progress = FALSE)
  
  # time column names differ; catch the usual ones
  nms <- names(df)
  tcol <- nms[which.max(grepl("time|timestamp|date", tolower(nms)))]
  if (is.na(tcol) || tcol == "") stop("GeoSphere CSV: no obvious time column. Columns: ", paste(nms, collapse=", "))
  
  # pick the TL column explicitly (avoid accidentally taking station_id=38)
  val_candidates <- setdiff(nms, tcol)
  val_col <- val_candidates[which(tolower(val_candidates) %in% c("tl", "temperature", "air_temperature"))[1]]
  if (is.na(val_col) || !nzchar(val_col)) {
    # common in GeoSphere CSV: parameter column is literally 'tl'
    # if not found, take the LAST non-time column (usually the parameter)
    val_col <- tail(val_candidates, 1)
  }
  
  out <- df |>
    transmute(
      timestamp = parse_dt_safe(.data[[tcol]], tz = "UTC") |> with_tz(TZ_LOCAL),
      TL = to_num(.data[[val_col]])
    ) |>
    filter(!is.na(timestamp)) |>
    arrange(timestamp) |>
    filter(timestamp >= as.POSIXct(start_date, tz = TZ_LOCAL),
           timestamp <  as.POSIXct(end_date + 1, tz = TZ_LOCAL))
  
  if (nrow(out) == 0) stop("GeoSphere TL: no rows after parsing/filter.")
  
  message("GeoSphere TL ok: station=", station_id,
          " rows=", nrow(out),
          " [", format(min(out$timestamp)), " .. ", format(max(out$timestamp)), "]")
  
  out
}

# ----------------------------
# LWD TL (LT): seasonal first, then *_latest
# - Reuse the proven logic from diagram_uid.R (httr2 + resp_body_raw)
# - Add fallback using *_LT_latest.csv
# ----------------------------
season_label <- function(date) {
  y <- year(date); m <- month(date)
  ifelse(m >= 10, sprintf("%d_%d", y, y + 1), sprintf("%d_%d", y - 1, y))
}

# Switch: Sys.setenv(DEBUG_LWD="1")
DEBUG_LWD <- identical(Sys.getenv("DEBUG_LWD", "0"), "1")

read_lwd_url <- function(url) {
  # Some clients occasionally receive 400/403 from wiski, so:
  # - accept content types
  # - set Referer
  # - retry
  resp <- tryCatch(
    request(url) |>
      req_headers(
        "Accept" = "text/csv,text/plain,*/*",
        "Referer" = "https://wiski.tirol.gv.at/",
        "Accept-Language" = "de,en;q=0.8"
      ) |>
      req_user_agent("icefall-model/1.0 (R httr2)") |>
      req_retry(max_tries = 4) |>
      req_error(is_error = function(r) FALSE) |>
      req_perform(),
    error = function(e) NULL
  )
  
  if (is.null(resp)) return(NULL)
  st <- resp_status(resp)
  
  raw <- tryCatch(resp_body_raw(resp), error = function(e) NULL)
  if (DEBUG_LWD) {
    msg <- tryCatch(resp_body_string(resp), error = function(e) "")
    msg <- substr(msg, 1, 200)
    message("LWD GET ", st, " -> ", url)
    if (nzchar(msg)) message("LWD body head: ", gsub("
", " ", msg))
  }
  
  if (st >= 400 || is.null(raw) || length(raw) == 0) return(NULL)
  
  txt0 <- tryCatch(rawToChar(raw), error = function(e) NA_character_)
  if (is.na(txt0)) return(NULL)
  txt <- iconv(txt0, from = "", to = "UTF-8", sub = "byte")
  
  # LWD files typically use ';' as delimiter
  tmp <- tryCatch(readr::read_delim(I(txt), delim = ";", show_col_types = FALSE, progress = FALSE),
                  error = function(e) NULL)
  if (is.null(tmp) || nrow(tmp) == 0) return(NULL)
  attr(tmp, "_source_url") <- url
  tmp
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
  df <- read_lwd_url(url)
  if (is.null(df)) return(NULL)
  
  nms <- names(df)
  nms_low <- tolower(iconv(nms, from = "", to = "ASCII//TRANSLIT", sub = ""))
  dt_i <- which(grepl("datetime|date_time|zeitstempel|timestamp|datumzeit|datum_zeit|zeit|datum", nms_low))
  dt_col <- if (length(dt_i)) nms[dt_i[1]] else NA_character_
  
  if (!is.na(dt_col)) {
    t <- parse_dt_any(df[[dt_col]], tz = TZ_LOCAL)
    val_col <- setdiff(nms, dt_col)[1]
    return(tibble(timestamp = t, value = to_num(df[[val_col]])) %>% filter(!is.na(timestamp)))
  }
  
  # fallback: first two cols
  if (ncol(df) >= 2) {
    t <- parse_dt_any(df[[1]], tz = TZ_LOCAL)
    return(tibble(timestamp = t, value = to_num(df[[2]])) %>% filter(!is.na(timestamp)))
  }
  NULL
}

read_lwd_latest <- function(station_code, param = "LT") {
  url <- sprintf("https://wiski.tirol.gv.at/lawine/produkte/ogd/%s/%s_%s_latest.csv",
                 station_code, station_code, param)
  df <- read_lwd_url(url)
  if (is.null(df)) return(NULL)
  
  nms <- names(df)
  nms_low <- tolower(iconv(nms, from = "", to = "ASCII//TRANSLIT", sub = ""))
  dt_i <- which(grepl("datetime|date_time|zeitstempel|timestamp|datumzeit|datum_zeit|zeit|datum", nms_low))
  dt_col <- if (length(dt_i)) nms[dt_i[1]] else NA_character_
  
  if (!is.na(dt_col)) {
    t <- parse_dt_any(df[[dt_col]], tz = TZ_LOCAL)
    val_col <- setdiff(nms, dt_col)[1]
    return(tibble(timestamp = t, value = to_num(df[[val_col]])) %>% filter(!is.na(timestamp)))
  }
  
  if (ncol(df) >= 2) {
    t <- parse_dt_any(df[[1]], tz = TZ_LOCAL)
    return(tibble(timestamp = t, value = to_num(df[[2]])) %>% filter(!is.na(timestamp)))
  }
  NULL
}

get_lwd_tl <- function(start_date, end_date, station_code) {
  # 1) seasonal (Oct 1 through today)
  seasons <- unique(season_label(seq(as.Date(start_date), as.Date(end_date), by = "day")))
  tl <- bind_rows(lapply(seasons, function(seas) read_lwd_param(station_code, "LT", seas)))
  
  # 2) fallback: latest (if the station lacks seasonal files)
  if (nrow(tl) == 0) tl <- read_lwd_latest(station_code, "LT")
  
  if (is.null(tl) || nrow(tl) == 0) {
    stop("No LWD TL data for station ", station_code,
         " (seasonal + latest failed). Tip: Sys.setenv(DEBUG_LWD=\"1\")")
  }
  
  out <- tl %>%
    transmute(timestamp, TL = value) %>%
    arrange(timestamp) %>%
    filter(timestamp >= as.POSIXct(start_date, tz = TZ_LOCAL),
           timestamp <  as.POSIXct(end_date + 1, tz = TZ_LOCAL))
  
  if (nrow(out) == 0) stop("LWD TL: 0 rows after date filter. station=", station_code)
  
  message("LWD TL ok: station=", station_code,
          " rows=", nrow(out),
          " [", format(min(out$timestamp)), " .. ", format(max(out$timestamp)), "]")
  
  out
}

# ----------------------------
# 1) Load
# ----------------------------
# IDs:
#  - GeoSphere: "38"  (Imst)
#  - LWD:      "IMUT2" (Muttekopfhütte)
#  - LWD:      "IMUT1" (Vorderer Maldonkopf)

imst <- get_geosphere_tl(START_DATE, END_DATE, station_id = "38") %>% rename(T0 = TL)  # 860m
mutk <- get_lwd_tl(START_DATE, END_DATE, station_code = "IMUT2")  %>% rename(T1 = TL)  # 1935m
mald <- get_lwd_tl(START_DATE, END_DATE, station_code = "IMUT1")  %>% rename(T2 = TL)  # 2580m

# ----------------------------
# 2) Align to 10-min grid + join
# ----------------------------
min_ts <- min(imst$timestamp, mutk$timestamp, mald$timestamp, na.rm = TRUE)
max_ts <- max(imst$timestamp, mutk$timestamp, mald$timestamp, na.rm = TRUE)

grid <- tibble(time = seq(
  from = floor_date(min_ts, step_str),
  to   = ceiling_date(max_ts, step_str),
  by   = step_str
))

roll_to_10min <- function(df, value_col, out_name) {
  df %>%
    mutate(time = floor_date(timestamp, step_str)) %>%
    group_by(time) %>%
    summarise(val = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    right_join(grid, by = "time") %>%
    arrange(time) %>%
    mutate(
      val = zoo::na.locf(val, na.rm = FALSE),
      val = zoo::na.locf(val, na.rm = FALSE, fromLast = TRUE)
    ) %>%
    transmute(time, !!out_name := val)
}

df0 <- roll_to_10min(imst, value_col = "T0", out_name = "T0")
df1 <- roll_to_10min(mutk, value_col = "T1", out_name = "T1")
df2 <- roll_to_10min(mald, value_col = "T2", out_name = "T2")

inv <- df0 %>%
  left_join(df1, by = "time") %>%
  left_join(df2, by = "time") %>%
  mutate(
    # Temperature differences (deg C)
    I01 = T1 - T0,
    I12 = T2 - T1,
    I02 = T2 - T0,
    
    # Physical gradients dT/dz (K/m); positive = inversion
    grad01_K_per_m = I01 / (Z1_m - Z0_m),
    grad12_K_per_m = I12 / (Z2_m - Z1_m),
    grad02_K_per_m = I02 / (Z2_m - Z0_m),
    
    # Express inversion strength as a positive gradient (°C/100 m)
    inv_grad_max_C_per_100m = pmax(0, pmax(grad01_K_per_m, grad12_K_per_m)) * 100,
    
    inv_score_C = 0.7 * pmax(0, I01) + 0.3 * pmax(0, I12),
    
    inv_class = case_when(
      inv_grad_max_C_per_100m < 0.05 ~ "none",
      inv_grad_max_C_per_100m < 0.20 ~ "weak",
      inv_grad_max_C_per_100m < 0.50 ~ "moderate",
      TRUE                            ~ "strong"
    )
  )

# Persistence/hysteresis
win_on  <- as.integer(3 * 60 / MODEL_STEP_MIN)  # 18 steps
win_off <- as.integer(2 * 60 / MODEL_STEP_MIN)  # 12 steps

roll_all <- function(x, w) zoo::rollapplyr(x, w, FUN = function(v) all(v, na.rm = TRUE), fill = FALSE, partial = TRUE)

cond_on  <- inv$inv_grad_max_C_per_100m >= 0.20
cond_off <- inv$inv_grad_max_C_per_100m <  0.05

# Persistence/hysteresis
win_on  <- as.integer(3 * 60 / MODEL_STEP_MIN)  # 18 steps
win_off <- as.integer(2 * 60 / MODEL_STEP_MIN)  # 12 steps

roll_all <- function(x, w) zoo::rollapplyr(x, w, FUN = function(v) all(v, na.rm = TRUE),
                                           fill = FALSE, partial = TRUE)

cond_on  <- inv$inv_grad_max_C_per_100m >= 0.20   # Inversion "an"
cond_off <- inv$inv_grad_max_C_per_100m <  0.05   # Inversion "aus"

on_ok  <- roll_all(cond_on,  win_on)
off_ok <- roll_all(cond_off, win_off)

inv_active <- logical(nrow(inv))
state <- FALSE
for (i in seq_len(nrow(inv))) {
  if (!state && on_ok[i])  state <- TRUE
  if ( state && off_ok[i]) state <- FALSE
  inv_active[i] <- state
}

inv <- inv %>%
  mutate(
    inv_active = inv_active,
    Z0_m = Z0_m, Z1_m = Z1_m, Z2_m = Z2_m
  ) %>%
  select(time, inv_active, inv_class,
         inv_grad_max_C_per_100m, inv_score_C,
         grad01_K_per_m, grad12_K_per_m, grad02_K_per_m,
         I01, I12, I02, Z0_m, Z1_m, Z2_m)

on_ok  <- roll_all(cond_on,  win_on)
off_ok <- roll_all(cond_off, win_off)

inv_active <- logical(nrow(inv))
state <- FALSE
for (i in seq_len(nrow(inv))) {
  if (!state && on_ok[i])  state <- TRUE
  if ( state && off_ok[i]) state <- FALSE
  inv_active[i] <- state
}

inv <- inv %>%
  mutate(inv_active = inv_active) %>%
  select(
    time, inv_active, inv_class,
    inv_grad_max_C_per_100m, inv_score_C,
    grad01_K_per_m, grad12_K_per_m, grad02_K_per_m,
    I01, I12, I02
  )

# ----------------------------
# 3) Save
# ----------------------------
dir.create("data/_cache_inversion", recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path("data/_cache_inversion", sprintf("inversion_%s.rds", format(END_DATE, "%Y%m%d")))
saveRDS(inv, out_rds)

message("Inversion cache written: ", out_rds)
