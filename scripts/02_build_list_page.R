# ============================================================
# Build list page (summary table) for GitHub Pages + offline viewing
# - meta:        data/Koordinaten_Wasserfaelle/eisklettern_links_entries_diff.csv
# - assignments: data/AWS/icefalls_nearest_station.csv (optional)
# - sun:         data/suntime/sun_uid_<uid>.csv (optional)
# - model runs:  data/ModelRuns/model_uid<uid>.csv
# - outputs:     site/icefalls_table.json + site/list.html
#   plus copies to repo root: icefalls_table.json + list.html
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(jsonlite)
  library(tibble)
})

TZ_LOCAL <- "Europe/Vienna"

# Current local day (used for sun summaries)
today <- as.Date(with_tz(Sys.time(), TZ_LOCAL))

# ----------------------------
# Paths
# ----------------------------
PATH_ASSIGN <- "data/AWS/icefalls_nearest_station.csv"
PATH_META   <- "data/Koordinaten_Wasserfaelle/eisklettern_links_entries_diff.csv"
DIR_SUN     <- "data/suntime"
DIR_MODELS  <- "data/ModelRuns"

OUT_DIR  <- "site"
OUT_JSON <- file.path(OUT_DIR, "icefalls_table.json")
OUT_HTML <- file.path(OUT_DIR, "list.html")

ROOT_JSON <- "icefalls_table.json"
ROOT_HTML <- "list.html"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helpers
# ----------------------------
parse_uid <- function(x) {
  as.integer(readr::parse_number(as.character(x)))
}

to_num <- function(x) {
  if (is.null(x)) return(NA_real_)
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

read_any_delim <- function(path, force_character = FALSE) {
  col_spec <- if (isTRUE(force_character)) readr::cols(.default = readr::col_character()) else readr::cols()
  x <- tryCatch(readr::read_delim(path, delim = "\t", col_types = col_spec, show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (!is.null(x) && ncol(x) > 1) return(x)
  x <- tryCatch(readr::read_delim(path, delim = ";", col_types = col_spec, show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (!is.null(x) && ncol(x) > 1) return(x)
  readr::read_csv(path, col_types = col_spec, show_col_types = FALSE, progress = FALSE)
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
    df <- tryCatch(read_any_delim(f), error = function(e) NULL)
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

get_chr <- function(df, ...) {
  cands <- c(...)
  for (nm in cands) if (nm %in% names(df)) return(as.character(df[[nm]]))
  rep(NA_character_, nrow(df))
}

get_num <- function(df, ...) {
  cands <- c(...)
  for (nm in cands) if (nm %in% names(df)) return(to_num(df[[nm]]))
  rep(NA_real_, nrow(df))
}

first_nonempty <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  x[[1]]
}

parse_time_any <- function(x, tz = TZ_LOCAL) {
  if (inherits(x, "POSIXct")) return(with_tz(x, tz))
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_

  # Parse timezone-aware timestamps (Z / +HH:MM) in UTC first, then convert.
  has_tz_offset <- grepl("(Z|[+-][0-9]{2}:?[0-9]{2})$", x)
  t <- as.POSIXct(rep(NA_real_, length(x)), origin = "1970-01-01", tz = tz)
  if (any(has_tz_offset, na.rm = TRUE)) {
    x_tz <- x[has_tz_offset]
    t_tz <- suppressWarnings(lubridate::parse_date_time(
      x_tz,
      orders = c("Y-m-dTH:M:SZ", "Y-m-dTH:MZ", "Y-m-dTH:M:S z", "Y-m-dTH:M z", "Y-m-d H:M:S z", "Y-m-d H:M z"),
      tz = "UTC"
    ))
    t[has_tz_offset] <- suppressWarnings(with_tz(t_tz, tz))
  }

  no_tz <- !has_tz_offset
  if (any(no_tz, na.rm = TRUE)) {
    x_local <- x[no_tz]
    x_local <- gsub("T", " ", x_local, fixed = TRUE)
    t_local <- suppressWarnings(lubridate::ymd_hms(x_local, tz = tz))
    if (all(is.na(t_local))) t_local <- suppressWarnings(lubridate::ymd_hm(x_local, tz = tz))
    if (all(is.na(t_local))) t_local <- suppressWarnings(lubridate::parse_date_time(
      x_local,
      orders = c("Ymd HMS", "Ymd HM", "Y-m-d H:M:S", "Y-m-d H:M", "Y-m-dTH:M:S", "Y-m-dTH:M"),
      tz = tz
    ))
    t[no_tz] <- t_local
  }

  t
}

parse_time_iso_z <- function(x, out_tz = TZ_LOCAL) {
  if (is.null(x)) return(as.POSIXct(NA))
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x2 <- sub("Z$", "", x)
  x2 <- gsub("T", " ", x2, fixed = TRUE)
  t <- suppressWarnings(lubridate::ymd_hms(x2, tz = "UTC"))
  if (all(is.na(t))) t <- suppressWarnings(lubridate::ymd_hm(x2, tz = "UTC"))
  suppressWarnings(lubridate::with_tz(t, out_tz))
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), NA_character_)
}

fmt_pct <- function(x, digits = 0) {
  ifelse(is.finite(x), paste0(round(x * 100, digits), "%"), NA_character_)
}

fmt_hm <- function(t) {
  ifelse(!is.na(t), format(t, "%H:%M"), NA_character_)
}

fmt_duration_h <- function(h) {
  h <- to_num(h)
  out <- rep(NA_character_, length(h))
  ok <- is.finite(h)
  mins <- round(h[ok] * 60)
  hh <- mins %/% 60
  mm <- mins %% 60
  out[ok] <- paste0(hh, " h ", mm, " min")
  out
}

normalize_angle_deg <- function(x) {
  ((x %% 360) + 360) %% 360
}

angle_diff_deg <- function(a, b) {
  a <- normalize_angle_deg(a)
  b <- normalize_angle_deg(b)
  d <- abs(a - b)
  pmin(d, 360 - d)
}

normalize_cardinal <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("NORTH", "N", x, fixed = TRUE)
  x <- gsub("SOUTH", "S", x, fixed = TRUE)
  x <- gsub("EAST", "E", x, fixed = TRUE)
  x <- gsub("WEST", "W", x, fixed = TRUE)
  x <- gsub("[^NSEW]", "", x)
  x[x == "NO"] <- "NE"
  x[x == "O"] <- "E"
  x[x == "SO"] <- "SE"
  x[x == "NNO"] <- "NNE"
  x[x == "ONO"] <- "ENE"
  x[x == "OSO"] <- "ESE"
  x[x == "SSO"] <- "SSE"
  x
}

cardinal_to_degree <- function(x) {
  key <- normalize_cardinal(x)
  lut <- c(
    "N" = 0, "NNE" = 22.5, "NE" = 45, "ENE" = 67.5,
    "E" = 90, "ESE" = 112.5, "SE" = 135, "SSE" = 157.5,
    "S" = 180, "SSW" = 202.5, "SW" = 225, "WSW" = 247.5,
    "W" = 270, "WNW" = 292.5, "NW" = 315, "NNW" = 337.5
  )
  out <- unname(lut[key])
  out[is.na(out)] <- NA_real_
  as.numeric(out)
}

is_degree_consistent_with_cardinal <- function(deg, cardinal) {
  target <- cardinal_to_degree(cardinal)
  is.finite(deg) & is.finite(target) & angle_diff_deg(deg, target) <= 45
}

# Fix common encoding / symbol issues (e.g., degree sign)
normalize_text <- function(x) {
  if (is.null(x)) return(NA_character_)
  x <- as.character(x)
  # best-effort to UTF-8 (no-op if already UTF-8)
  x <- suppressWarnings(iconv(x, from = "", to = "UTF-8"))
  # normalize degree symbol variants / mojibake
  x <- gsub("Â°", "°", x, fixed = TRUE)
  x <- gsub("º",  "°", x, fixed = TRUE)
  x <- gsub("ø",  "°", x, fixed = TRUE)
  x
}



# ----------------------------
# 1) Load meta (CSV)
# ----------------------------
if (!file.exists(PATH_META)) stop("Missing: ", PATH_META)

meta_raw <- read_any_delim(PATH_META, force_character = TRUE) %>%
  rename_with(tolower)

if (!"uid" %in% names(meta_raw)) stop("META CSV has no 'uid' column.")

meta <- tibble(
  uid = parse_uid(meta_raw$uid),
  name = get_chr(meta_raw, "name"),
  topo_url = get_chr(meta_raw, "topo_url"),
  topo_slug = get_chr(meta_raw, "topo_slug"),
  latitude  = get_num(meta_raw, "latitude", "lat"),
  longitude = get_num(meta_raw, "longitude", "lon"),
  elev_m = get_num(meta_raw, "elevation_dgm5m", "elevation_m", "elevation", "elev_m", "altitude_m"),
  difficulty = get_chr(meta_raw, "difficulty", "grade"),
  aspect = get_num(meta_raw, "aspect_deg", "aspect", "exposure_deg"),
  aspect_deg = get_num(meta_raw, "aspect_deg", "aspect", "exposure_deg"),
  aspect_cardinal = normalize_cardinal(get_chr(meta_raw, "aspect_cardinal", "exposure", "direction")),
  icefall_height_m = get_num(meta_raw, "icefall_height_m", "height_m", "height"),
  approach = get_chr(meta_raw, "approach"),
  descent  = get_chr(meta_raw, "descent"),
  first_ascent = get_chr(meta_raw, "first_ascent"),
  description  = get_chr(meta_raw, "description")
) %>%
  filter(!is.na(uid)) %>%
  mutate(dplyr::across(where(is.character), normalize_text))

uid_raw <- as.character(meta_raw$uid)
uid_parsed <- parse_uid(uid_raw)
message(
  "META: rows=", nrow(meta_raw),
  " | parsed rows=", nrow(meta),
  " | unique uids=", length(unique(meta$uid)),
  " | uid NA=", sum(is.na(uid_parsed)), "/", length(uid_parsed)
)
message("META uid examples: ", paste(head(uid_raw, 8), collapse = " | "))
if (nrow(meta) <= 1) {
  warning(
    "META parsed to ", nrow(meta), " row(s).",
    " This often means PATH_META is a diff/partial file or uid parsing failed."
  )
}

# ----------------------------
# 2) Load assign (optional)
# ----------------------------
assign <- NULL
if (file.exists(PATH_ASSIGN)) {
  assign <- read_any_delim(PATH_ASSIGN) %>%
    mutate(uid = parse_uid(uid))
}

# ----------------------------
# 3) Sun horizons (optional)
# ----------------------------
sun <- NULL
sun_missing_uids <- integer(0)
if (dir.exists(DIR_SUN)) {
  sun_loaded <- load_sun_for_uids(unique(meta$uid), dir_sun = DIR_SUN)
  sun_missing_uids <- sun_loaded$missing_uids
  if (length(sun_loaded$missing_uids) > 0) {
    message(
      "Warning: missing sun files for UIDs: ",
      paste(sprintf("%03d", sun_loaded$missing_uids), collapse = ", ")
    )
  }

  sun_raw <- sun_loaded$data %>%
    rename_with(tolower) %>%
    mutate(
      uid = parse_uid(uid),
      date = as.Date(get_chr(., "date"))
    )

  if (!"sun_hours_topo" %in% names(sun_raw)) {
    sun_raw$sun_hours_topo <- NA_real_
  }
  
  sun <- sun_raw %>%
    filter(.data$date == today) %>%
    group_by(uid) %>%
    summarise(
      sunrise_topo_local = parse_time_iso_z(first_nonempty(sunrise_topo), out_tz = TZ_LOCAL),
      sunset_topo_local  = parse_time_iso_z(first_nonempty(sunset_topo),  out_tz = TZ_LOCAL),
      sun_hours_today_h = to_num(first_nonempty(sun_hours_topo)),
      sun_today_range_txt = dplyr::if_else(
        is.na(sunrise_topo_local) | is.na(sunset_topo_local),
        NA_character_,
        paste0(fmt_hm(sunrise_topo_local), "-", fmt_hm(sunset_topo_local))
      ),
      sun_duration_today_txt = fmt_duration_h(sun_hours_today_h),
      .groups = "drop"
    )
  
  if (all(is.na(sun$sun_hours_today_h))) {
    sun <- sun %>%
      mutate(
        sun_hours_today_h = as.numeric(difftime(sunset_topo_local, sunrise_topo_local, units = "hours")),
        sun_duration_today_txt = fmt_duration_h(sun_hours_today_h)
      )
  }
} else {
  message("Warning: sun directory missing: ", DIR_SUN)
}

# ----------------------------
# 4) Model summary (latest historical value)
# ----------------------------
summarise_uid_model <- function(uid) {
  f <- file.path(DIR_MODELS, sprintf("model_uid%s.csv", uid))
  
  empty <- tibble(
    uid = uid,
    thickness_latest_m = NA_real_,
    thickness_latest_time_local = NA_character_
  )
  
  if (!file.exists(f) || file.info(f)$size <= 0) return(empty)
  
  df <- tryCatch(readr::read_csv(f, show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (is.null(df) || !"time" %in% names(df)) return(empty)
  
  df <- df %>%
    mutate(
      time = parse_time_any(.data$time, tz = TZ_LOCAL),
      date = as.Date(with_tz(time, TZ_LOCAL)),
      thickness_m = if ("thickness_m" %in% names(df)) to_num(.data$thickness_m) else NA_real_
    ) %>%
    filter(!is.na(time), is.finite(thickness_m)) %>%
    arrange(time)
  
  if (nrow(df) == 0) return(empty)
  latest <- df[nrow(df), , drop = FALSE]
  
  tibble(
    uid = uid,
    thickness_latest_m = latest$thickness_m[[1]],
    thickness_latest_time_local = format(with_tz(latest$time[[1]], TZ_LOCAL), "%d.%m.%Y %H:%M")
  )
}

uids <- sort(unique(meta$uid))
model_sum <- bind_rows(lapply(uids, summarise_uid_model))

# ----------------------------
# 5) Merge
# ----------------------------
out <- meta %>% left_join(model_sum, by = "uid")

if (!is.null(assign)) {
  assign_slim <- assign %>%
    dplyr::select(dplyr::any_of(c(
      "uid", "station_id", "source", "dist_km", "elev_diff_m",
      "icefall_name", "ice_lon", "ice_lat", "icefall_elev_m", "icefall_height_m"
    ))) %>%
    dplyr::mutate(uid = parse_uid(uid))
  
  out <- out %>%
    left_join(assign_slim, by = "uid") %>%
    mutate(
      name = coalesce(as.character(.data$icefall_name), .data$name, paste0("UID ", .data$uid)),
      latitude  = coalesce(to_num(.data$ice_lat), .data$latitude),
      longitude = coalesce(to_num(.data$ice_lon), .data$longitude),
      elev_m = coalesce(to_num(.data$icefall_elev_m), .data$elev_m),
      icefall_height_m = coalesce(to_num(.data$icefall_height_m), .data$icefall_height_m)
    )
} else {
  out$station_id <- NA_character_
  out$source <- NA_character_
  out$dist_km <- NA_real_
  out$elev_diff_m <- NA_real_
}

if (!is.null(sun)) {
  out <- out %>% left_join(sun, by = "uid")
} else {
  if (!"topo_url" %in% names(out)) out$topo_url <- NA_character_
  if (!"topo_slug" %in% names(out)) out$topo_slug <- NA_character_
  out$sun_today_range_txt <- NA_character_
  out$sun_hours_today_h <- NA_real_
  out$sun_duration_today_txt <- NA_character_
}

out <- out %>%
  mutate(
    sun_status = dplyr::case_when(
      uid %in% sun_missing_uids ~ "missing_file",
      is.finite(sun_hours_today_h) & !is.na(sun_today_range_txt) ~ "has_values",
      TRUE ~ "no_values_day"
    ),
    sun_today_range_txt = dplyr::case_when(
      sun_status == "missing_file" ~ "NO SUN DATA",
      sun_status == "no_values_day" ~ "no direct sunlight",
      TRUE ~ sun_today_range_txt
    )
  )

if ("topo_url.y" %in% names(out) || "topo_slug.y" %in% names(out)) {
  out <- out %>%
    mutate(
      topo_url = dplyr::coalesce(.data$topo_url, .data$topo_url.y),
      topo_slug = dplyr::coalesce(.data$topo_slug, .data$topo_slug.y)
    ) %>%
    select(-dplyr::any_of(c("topo_url.y", "topo_slug.y")))
}

out <- out %>%
  mutate(
    plot_url = sprintf("plots/uid_%03d.png", uid),
    thickness_latest_txt = fmt_num(thickness_latest_m, 2)
  ) %>%
  arrange(desc(thickness_latest_m))

# ----------------------------
# 6) Write JSON
# ----------------------------
jsonlite::write_json(out, OUT_JSON, pretty = TRUE, auto_unbox = TRUE, na = "null")
message("✅ Wrote JSON: ", OUT_JSON)

# ----------------------------
# 7) Write list.html
#    JSON is embedded as Base64, so the page also works offline (file://)
# ----------------------------
today_str <- format(today, "%d.%m.%Y")

embedded_json <- jsonlite::toJSON(out, auto_unbox = TRUE, na = "null")
embedded_b64  <- jsonlite::base64_enc(charToRaw(enc2utf8(embedded_json)))

# External JS (kept separate for maintainability)
ASSET_DIR <- file.path(OUT_DIR, "assets")
OUT_JS <- file.path(ASSET_DIR, "list.js")
dir.create(ASSET_DIR, showWarnings = FALSE, recursive = TRUE)

js_source_path <- file.path('scripts', 'list.js')
js_lines <- readLines(js_source_path, warn = FALSE)
if (!length(js_lines)) {
  stop('Missing JS source: ', js_source_path)
}

# Build HTML as lines to avoid locale/size parser limits (notably on Windows)
# and to keep the R source ASCII-only (umlauts via HTML entities).
html_lines <- c(
  '<!doctype html>',
  '<html lang="en">',
  '<head>',
  '  <meta charset="utf-8"/>',
  '  <meta name="viewport" content="width=device-width, initial-scale=1"/>',
  '  <title>Icefalls - Overview</title>',
  '  <style>',
  '    body { font-family: system-ui, -apple-system, Segoe UI, Roboto, sans-serif; margin: 0; padding: 0; color: #111; }',
  '    header { padding: 10px 14px; border-bottom: 1px solid #ddd; display:flex; gap:12px; align-items:center; flex-wrap:wrap; }',
  '    header a { text-decoration:none; padding:6px 10px; border:1px solid #ddd; border-radius:8px; color:#111; }',
  '    header a:hover { background:#f4f4f4; }',
  '    .wrap { padding: 12px 14px; }',
  '    .controls { display:flex; gap:10px; flex-wrap:wrap; align-items:center; margin-bottom:10px; }',
  '    .controls > * { flex: 0 1 auto; }',
  '    input[type="search"], input[type="number"], select { padding:10px 12px; border:1px solid #ccc; border-radius:10px; font-size:16px; }',
  '    input[type="range"] { width: 220px; }',
  '    .table-wrap { overflow-x:auto; border:1px solid #eee; border-radius:12px; }',
  '    table { width:100%; min-width: 980px; border-collapse: collapse; }',
  '    th, td { padding: 10px 8px; border-bottom: 1px solid #eee; vertical-align: top; }',
  '    th { text-align:left; position: sticky; top: 0; background: #fff; z-index: 1; cursor:pointer; user-select:none; }',
  '    tr:hover { background: #fafafa; }',
  '    .muted { color:#666; font-size:12px; }',
  '    .btn { display:inline-flex; gap:6px; align-items:center; padding:6px 10px; border:1px solid #ddd; border-radius:10px; background:#fff; cursor:pointer; }',
  '    .btn:hover { background:#f4f4f4; }',
  '    .small { font-size: 12px; }',
  '    details > summary { list-style: none; }',
  '    details > summary::-webkit-details-marker { display:none; }',
  '    #filters .panel { margin-top:10px; padding:12px; border:1px solid #eee; border-radius:12px; background:#fafafa; display:flex; flex-direction:column; gap:12px; }',
  '    #filters .section { padding:10px; border:1px solid #e6e6e6; border-radius:12px; background:#fff; display:flex; flex-direction:column; gap:8px; }',
  '    #filters .section-title { font-size:12px; color:#444; font-weight:700; letter-spacing:0.02em; text-transform:uppercase; }',
  '    #filters .row { display:flex; flex-wrap:wrap; gap:10px; align-items:center; }',
  '    #filters label { font-size:12px; color:#666; min-width:140px; }',
  '    #filters input[type="text"] { min-width:240px; }',
  '    #modal { display:none; position:fixed; inset:0; background:rgba(0,0,0,0.8); z-index:9999; }',
  '    #modal .inner { position:absolute; inset:0; display:flex; flex-direction:column; }',
  '    #modal .bar { padding:10px; display:flex; gap:10px; align-items:center; justify-content:space-between; color:#fff; }',
  '    #modal img { flex:1; width:100%; height:100%; object-fit: contain; }',
  '    #modal .bar button, #modal .bar a {',
  '      color:#fff; border:1px solid rgba(255,255,255,0.35);',
  '      background: transparent; padding:8px 12px; border-radius:10px; cursor:pointer;',
  '      text-decoration:none;',
  '    }',
  '    #modal .bar button:hover, #modal .bar a:hover { background: rgba(255,255,255,0.12); }',
  '    @media (max-width: 720px) {',
  '      header { gap:8px; }',
  '      header a { flex: 1 1 auto; text-align:center; }',
  '      .controls { flex-direction: column; align-items: stretch; }',
  '      .controls > * { width: 100%; }',
  '      #filters .panel { padding:10px; }',
  '      #filters .section { gap:10px; }',
  '      #filters .row { flex-direction: column; align-items: stretch; }',
  '      #filters label { min-width: 0; }',
  '      input[type="search"] { width: 100%; min-width: 0; }',
  '      input[type="range"] { width: 100%; }',
  '      th, td { padding: 12px 6px; }',
  '      table { min-width: 760px; }',
  '    }',
  '  </style>',
  '</head>',
  '<body>',
  '  <header>',
  '    <a href="index.html">🏠 Home</a>',
  '    <a href="map.html">Map</a>',
  '    <a href="list.html"><b>Overview</b></a>',
  paste0('    <span class="muted">Today: ', today_str, ' (TZ: Europe/Vienna)</span>'),
  '  </header>',
  '',
  '  <div class="wrap">',
  '    <div class="controls">',
  '      <input id="q" type="search" placeholder="Search: name, grade, station ...">',
  '',
  '      <details id="filters">',
  '        <summary class="btn" type="button">Filters</summary>',
  '        <div class="panel">',
  '          <div class="section">',
  '            <div class="section-title">Location</div>',
  '            <div class="row">',
  '              <label for="radiusKm">Radius</label>',
  '              <input id="radiusKm" type="number" min="0" step="1" value="0" style="width:110px;" title="Radius in km (0 = off)">',
  '              <span class="muted">0 km = no filter</span>',
  '            </div>',
  '',
  '            <div class="row">',
  '              <label>Center</label>',
  '              <button class="btn" id="useGeo" type="button">GPS</button>',
  '              <input id="place" type="text" list="placeSuggestions" placeholder="Enter place (e.g. Obergurgl)">',
  '              <datalist id="placeSuggestions"></datalist>',
  '              <button class="btn" id="geocodeBtn" type="button">Search</button>',
  '            </div>',
  '',
  '            <div class="row">',
  '              <label>Coord.</label>',
  '              <input id="centerLat" type="number" step="0.000001" placeholder="lat" style="width:140px;">',
  '              <input id="centerLon" type="number" step="0.000001" placeholder="lon" style="width:140px;">',
  '              <button class="btn" id="setCustom" type="button">Set</button>',
  '            </div>',
  '',
  '            <div class="muted small">Note: GPS usually works only via https/localhost. Place search needs internet.</div>',
  '            <div class="muted small" id="geoStatus"></div>',
  '          </div>',
  '',
  '          <div class="section">',
  '            <div class="section-title" >Difficulty</div>',
  '            <div class="row">',
  '              <label title="Technical climbing (A) – e.g. A1 to A4">A</label>',
  '              <input id="aMin" type="range" min="0.75" max="4.25" step="0.25" value="0.75" title="A1- to A4+">',
  '              <input id="aMax" type="range" min="0.75" max="4.25" step="0.25" value="4.25" title="A1- to A4+">',
  '              <span class="muted" id="aRangeTxt">A1- – A4+</span>',
  '            </div>',
  '            <div class="row">',
  '              <label title="Mixed climbing (M) – combination of ice & rock">M</label>',
  '              <input id="mMin" type="range" min="0.75" max="13.25" step="0.25" value="0.75" title="M1- to M13+">',
  '              <input id="mMax" type="range" min="0.75" max="13.25" step="0.25" value="13.25" title="M1- to M13+">',
  '              <span class="muted" id="mRangeTxt">M1- – M13+</span>',
  '            </div>',
  '            <div class="row">',
  '              <label title="Water ice (WI) – pure ice grading">WI</label>',
  '              <input id="wiMin" type="range" min="0.75" max="7.25" step="0.25" value="0.75" title="WI1- to WI7+">',
  '              <input id="wiMax" type="range" min="0.75" max="7.25" step="0.25" value="7.25" title="WI1- to WI7+">',
  '              <span class="muted" id="wiRangeTxt">WI1- – WI7+</span>',
  '            </div>',
  '            <div class="row">',
  '              <label>Rock (UIAA)</label>',
  '              <input id="rMin" type="range" min="0.75" max="12.25" step="0.25" value="0.75" title="1- to 12+">',
  '              <input id="rMax" type="range" min="0.75" max="12.25" step="0.25" value="12.25" title="1- to 12+">',
  '              <span class="muted" id="rRangeTxt">1- – 12+</span>',
  '            </div>',
  '          </div>',
  '',
  '          <div class="section">',
  '            <div class="section-title">Sun</div>',
  '            <div class="row">',
  '              <label>Sun today (h)</label>',
  '              <input id="sunMin" type="range" min="0" max="12" step="0.25" value="0" title="Sun duration (topography)">',
  '              <input id="sunMax" type="range" min="0" max="12" step="0.25" value="12" title="Sun duration (topography)">',
  '              <span class="muted" id="sunRangeTxt">0.0 – 12.0 h</span>',
  '            </div>',
  '            <div class="row">',
  '              <label>Elevation (m)</label>',
  '              <input id="elevMin" type="range" min="0" max="4000" step="50" value="0" title="Icefall elevation in m">',
  '              <input id="elevMax" type="range" min="0" max="4000" step="50" value="4000" title="Icefall elevation in m">',
  '              <span class="muted" id="elevRangeTxt">0 – 4000 m</span>',
  '            </div>',
  '          </div>',
  '        </div>',
  '      </details>',
  '',
  '      <span class="muted" style="margin-left:auto;">Click column header = sort</span>',
  '    </div>',
  '',
  '    <div class="muted small" id="status">Loading data ...</div>',
  '',
  '    <div class="table-wrap">',
  '      <table id="tbl">',
  '        <thead>',
  '          <tr>',
  '            <th data-key="name">Icefall</th>',
  '            <th data-key="difficulty">Difficulty</th>',
  '            <th data-key="_grade_a" title="A = aid climbing">A</th>',
  '            <th data-key="_grade_m" title="M = mixed climbing (ice & rock)">M</th>',
  '            <th data-key="_grade_wi" title="WI = water ice">WI</th>',
  '            <th data-key="_grade_r">Rock</th>',
  '            <th data-key="elev_m">Elevation (m)</th>',
  '            <th data-key="_dist_km">Distance (km)</th>',
  '            <th data-key="sun_today_range_txt">Sun today</th>',
  '            <th data-key="sun_hours_today_h">Sun duration</th>',
  '            <th data-key="thickness_latest_m">Latest ice thickness (m)</th>',
  '            <th data-key="_last_upload_ts">Latest upload</th>',
  '            <th>Details</th>',
  '            <th>Topo</th>',
  '          </tr>',
  '        </thead>',
  '        <tbody></tbody>',
  '      </table>',
  '    </div>',
  '  </div>',
  '',
  '  <div id="modal">',
  '    <div class="inner">',
  '      <div class="bar">',
  '        <div id="modalTitle">Chart</div>',
  '        <div style="display:flex; gap:10px; align-items:center;">',
  '          <a id="openNewTab" href="#" target="_blank" rel="noopener">In new tab</a>',
  '          <button id="closeModal">Close</button>',
  '        </div>',
  '      </div>',
  '      <img id="modalImg" src="" alt="Chart"/>',
  '    </div>',
  '  </div>',
  '',
  paste0('  <script id="ICEFALL_DATA_B64" type="text/plain">', embedded_b64, '</script>'),
  '',
  '  <script src="assets/list.js"></script>',
  '</body>',
  '</html>'
)

html <- paste(html_lines, collapse = "\n")

writeLines(enc2utf8(js_lines), OUT_JS, useBytes = TRUE)
message("✅ Wrote JS: ", OUT_JS)

writeLines(enc2utf8(html), OUT_HTML, useBytes = TRUE)
message("✅ Wrote HTML: ", OUT_HTML)

# ----------------------------
# 8) Write detail pages per UID (site/icefalls/uid_###.html)
# ----------------------------

DETAIL_DIR <- file.path(OUT_DIR, "icefalls")
dir.create(DETAIL_DIR, showWarnings = FALSE, recursive = TRUE)
old_detail_pages <- list.files(DETAIL_DIR, pattern = "^uid_[0-9]+[.]html$", full.names = TRUE)
if (length(old_detail_pages) > 0) unlink(old_detail_pages)

# Cloudflare Worker + R2 Public base
API_BASE  <- "https://icefalls-api.carlos-wydra.workers.dev"
R2_PUBLIC <- "https://pub-1b553d6f009540c0881b434b7791c67c.r2.dev"

# Helper: safe HTML text (very small; JS does most escaping)
esc_html <- function(x) {
  x <- ifelse(is.na(x) | x %in% c("", "NA", "NaN", "NULL"), "", as.character(x))
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

# Iterate over all rows (one detail page per uid)
for (i in seq_len(nrow(out))) {
  
  r <- out[i, , drop = FALSE]
  
  uid <- as.integer(r$uid[[1]])
  if (!is.finite(uid)) next
  
  uid_pad <- sprintf("%03d", uid)
  
  nm  <- esc_html(r$name[[1]])
  diff <- esc_html(r$difficulty[[1]])

  aspect_cardinal_raw <- NA_character_
  if ("aspect_cardinal" %in% names(r) && length(r[["aspect_cardinal"]]) > 0) {
    aspect_cardinal_raw <- as.character(r[["aspect_cardinal"]][[1]])
  }
  aspect_cardinal_txt <- esc_html(aspect_cardinal_raw)

  aspect_deg <- NA_real_
  if ("aspect_deg" %in% names(r) && length(r[["aspect_deg"]]) > 0 && !is.na(r[["aspect_deg"]][[1]])) {
    aspect_deg <- suppressWarnings(as.numeric(r[["aspect_deg"]][[1]]))
  }
  if (!is.finite(aspect_deg)) {
    if ("aspect" %in% names(r) && length(r[["aspect"]]) > 0 && !is.na(r[["aspect"]][[1]])) {
      aspect_deg <- suppressWarnings(as.numeric(r[["aspect"]][[1]]))
      if (!is.finite(aspect_deg)) aspect_deg <- NA_real_
    } else {
      aspect_deg <- NA_real_
    }
  }

  # If source degree does not fit the cardinal direction, derive degree from cardinal.
  if (!is.na(aspect_cardinal_raw) &&
      nchar(trimws(aspect_cardinal_raw)) > 0 &&
      (is.na(aspect_deg) || !is.finite(aspect_deg) ||
       !is_degree_consistent_with_cardinal(aspect_deg, aspect_cardinal_raw))) {
    aspect_deg <- cardinal_to_degree(aspect_cardinal_raw)
  }
  aspect_deg_txt <- if (is.finite(aspect_deg)) paste0(round(aspect_deg), "&deg;") else ""
  aspect_txt <- dplyr::case_when(
    nchar(aspect_cardinal_txt) > 0 & nchar(aspect_deg_txt) > 0 ~ paste0(aspect_cardinal_txt, " (", aspect_deg_txt, ")"),
    nchar(aspect_cardinal_txt) > 0 ~ aspect_cardinal_txt,
    nchar(aspect_deg_txt) > 0 ~ aspect_deg_txt,
    TRUE ~ "&mdash;"
  )
  elev <- r$elev_m[[1]]
  elev_txt <- if (is.finite(elev)) paste0(round(elev), " m") else "&mdash;"
  
  topo_url <- as.character(r$topo_url[[1]])
  topo_url <- ifelse(is.na(topo_url) | topo_url == "", "", topo_url)
  
  lat <- r$latitude[[1]]
  lon <- r$longitude[[1]]
  lat_txt <- if (is.finite(lat)) sprintf("%.6f", lat) else ""
  lon_txt <- if (is.finite(lon)) sprintf("%.6f", lon) else ""

  sun_range_txt <- as.character(r$sun_today_range_txt[[1]])
  if (is.na(sun_range_txt) || sun_range_txt == "") sun_range_txt <- "&mdash;"
  sun_duration_txt <- as.character(r$sun_duration_today_txt[[1]])
  if (is.na(sun_duration_txt) || sun_duration_txt == "") sun_duration_txt <- "&mdash;"

  sun_series_json <- "[]"
  sun_file <- find_sun_file_for_uid(uid)
  if (!is.na(sun_file)) {
    sun_df_uid <- tryCatch(read_any_delim(sun_file), error = function(e) NULL)
    if (!is.null(sun_df_uid) && nrow(sun_df_uid) > 0) {
      sun_df_uid <- sun_df_uid %>%
        rename_with(tolower)
      if (all(c("date", "sun_hours_topo") %in% names(sun_df_uid))) {
        sun_df_uid <- sun_df_uid %>%
          transmute(
            date = as.Date(as.character(date)),
            sun_hours = to_num(sun_hours_topo)
          ) %>%
          filter(!is.na(date), is.finite(sun_hours)) %>%
          arrange(date)
        if (nrow(sun_df_uid) > 0) {
          sun_df_uid$date <- format(sun_df_uid$date, "%Y-%m-%d")
          sun_series_json <- jsonlite::toJSON(
            sun_df_uid,
            dataframe = "rows",
            auto_unbox = TRUE,
            na = "null"
          )
        }
      }
    }
  }
  
  # plot png relative to site root
  plot_rel <- sprintf("plots/uid_%03d.png", uid)
  
  out_file <- file.path(DETAIL_DIR, sprintf("uid_%s.html", uid_pad))
  
  detail_html <- c(
    "<!doctype html>",
    "<html lang='en'>",
    "<head>",
    "  <meta charset='utf-8'/>",
    "  <meta name='viewport' content='width=device-width, initial-scale=1'/>",
    paste0("  <title>", nm, " (UID ", uid_pad, ")</title>"),
    
    # Leaflet (small map)
    "  <link rel='stylesheet' href='https://unpkg.com/leaflet@1.9.4/dist/leaflet.css'/>",
    "  <script src='https://unpkg.com/leaflet@1.9.4/dist/leaflet.js'></script>",
    
    "  <style>",
    "    body{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;margin:0;padding:0;background:#fff;}",
    "    header{padding:10px 14px;border-bottom:1px solid #ddd;display:flex;gap:10px;align-items:center;flex-wrap:wrap;}",
    "    header a{text-decoration:none;padding:6px 10px;border:1px solid #ddd;border-radius:10px;color:#111;background:#fff;}",
    "    header a:hover{background:#f4f4f4;}",
    "    .wrap{padding:12px 14px;max-width:1400px;margin:0 auto;}",
    
    "    .grid2{display:grid;grid-template-columns:1fr 1fr;gap:12px;}",
    "    @media(max-width:900px){.grid2{grid-template-columns:1fr;}}",
    
    "    .card{border:1px solid #eee;border-radius:16px;padding:12px 14px;background:#fff;}",
    "    .card h2{margin:0 0 8px 0;font-size:16px;}",
    "    .kv{display:grid;grid-template-columns:170px 1fr;gap:6px 10px;font-size:14px;}",
    "    .k{color:#666;}",
    "    .muted{color:#666;font-size:12px;}",
    
    "    .btn{display:inline-flex;gap:6px;align-items:center;padding:8px 12px;border:1px solid #ddd;border-radius:12px;background:#fff;cursor:pointer;text-decoration:none;color:#111;}",
    "    .btn:hover{background:#f4f4f4;}",
    
    "    .gallery{display:flex;flex-wrap:wrap;gap:10px;margin-top:6px;}",
    "    .ph{width:220px;border:1px solid #eee;border-radius:14px;overflow:hidden;background:#fff;}",
    "    .ph img{width:100%;height:160px;object-fit:cover;display:block;}",
    "    .cap{padding:8px 10px;font-size:12px;}",
    
    "    input[type='file'], input[type='date'], select, textarea{padding:10px 12px;border:1px solid #ccc;border-radius:12px;font-size:14px;}",
    "    textarea{resize:vertical;}",
    
    # diagram full width
    "    .plotWrap{margin-top:12px;}",
    "    .plotImg{width:100%;height:auto;display:block;border:1px solid #eee;border-radius:16px;cursor:zoom-in;}",
    "    .chartWrap{margin-top:12px;}",
    "    .chartCanvas{width:100%;height:280px;display:block;border:1px solid #eee;border-radius:16px;background:#fff;}",
    "    .chartHint{margin-top:8px;font-size:12px;color:#666;}",
    
    "    #map{height:240px;border-radius:16px;border:1px solid #eee;}",
    "    @media(max-width:720px){",
    "      header a{flex:1 1 auto;text-align:center;}",
    "      .wrap{padding:10px 12px;}",
    "      .kv{grid-template-columns:1fr;}",
    "      .gallery{flex-direction:column;}",
    "      .ph{width:100%;}",
    "      .ph img{height:auto;}",
    "      #map{height:200px;}",
    "    }",
    "  </style>",
    "</head>",
    "<body>",
    
    "<header>",
    "  <a href='../index.html'>🏠 Home</a>",
    "  <a href='../list.html'>&larr; Overview</a>",
    "  <a href='../map.html'>Map</a>",
    paste0("  <div style='font-weight:700;'>", nm, " <span class='muted'>(UID ", uid_pad, ")</span></div>"),
    "</header>",
    
    "<div class='wrap'>",
    
    # top row: basic info + images (empty state initially)
    "  <div class='grid2'>",
    "    <div class='card'>",
    "      <h2>Basic info</h2>",
    "      <div class='kv'>",
    paste0("        <div class='k' >Difficulty</div><div><b>", ifelse(nchar(diff)>0, diff, "&mdash;"), "</b></div>"),
    paste0("        <div class='k'>Aspect</div><div><b>", aspect_txt, "</b></div>"),
    paste0("        <div class='k'>Elevation a.s.l.</div><div><b>", elev_txt, "</b></div>"),
    paste0("        <div class='k'>Sunlight today</div><div><b>", sun_range_txt, "</b> <span class='muted'>(", sun_duration_txt, ")</span></div>"),
    paste0("        <div class='k'>Topo</div><div>", ifelse(nchar(topo_url)>0, paste0("<a href='", topo_url, "' target='_blank' rel='noopener'>Topo</a>"), "<span class='muted'>&mdash;</span>"), "</div>"),
    "      </div>",
    "    </div>",
    
    "    <div class='card'>",
    "      <h2>Photos</h2>",
    "      <div class='muted' style='display:flex;justify-content:space-between;gap:10px;align-items:center;'>",
    "        <div id='imgStatus'>Loading...</div>",
    "        <div>Uploads are displayed publicly.</div>",
    "      </div>",
    "      <div id='gallery' class='gallery' style='margin-top:10px;'></div>",
    "    </div>",
    "  </div>",
    
    # full-width plot card
    "  <div class='card plotWrap'>",
    "    <h2>Chart</h2>",
    paste0("    <a href='../", plot_rel, "' target='_blank' rel='noopener' title='Open chart full size'>"),
    paste0("      <img class='plotImg' src='../", plot_rel, "' alt='Chart UID ", uid_pad, "' onerror=\"this.outerHTML='<div class=&quot;muted&quot;>No plot found.</div>';\"/>"),
    "    </a>",
    "    <div class='chartWrap'>",
    "      <h2>Seasonal sunlight</h2>",
    "      <canvas id='sunSeasonChart' class='chartCanvas' width='1200' height='280'></canvas>",
    "      <div id='sunChartHover' class='chartHint'>Hover to show the date.</div>",
    "    </div>",
    "  </div>",
    
    # upload + map row
    "  <div class='grid2' style='margin-top:12px;'>",
    
    # upload card
    "    <div class='card'>",
    "      <h2>Upload photos</h2>",
    "      <div class='muted' style='margin-bottom:10px;'>By uploading you confirm that you own the rights to the image. Please do not include personal data.</div>",
    
    "      <div style='display:flex;gap:10px;flex-wrap:wrap;align-items:center;'>",
    "        <div class='muted' style='min-width:110px;'>Ice condition</div>",
    "        <select id='rating'>",
    "          <option value='3'>good</option>",
    "          <option value='2'>marginal</option>",
    "          <option value='1'>poor</option>",
    "        </select>",
    "        <div class='muted' style='min-width:50px;'>Date</div>",
    "        <input id='shotDate' type='date'/>",
    "      </div>",
    
    # comment ONLY here (NOT in basic info)
    "      <div style='margin-top:10px;'>",
    "        <textarea id='comment' rows='3' placeholder='Comment (optional): ice state, conditions, hazards…' style='width:100%;max-width:100%;'></textarea>",
    "      </div>",
    
    "      <div style='display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin-top:10px;'>",
    "        <input id='file' type='file' accept='image/*' multiple style='position:absolute;left:-9999px;width:1px;height:1px;overflow:hidden;'/>",
    "        <label for='file' class='btn' style='cursor:pointer;'>Choose files</label>",
    "        <span id='fileStatus' class='muted'>No files selected</span>",
    "        <button id='btnUpload' class='btn' type='button'>Upload</button>",
    "        <span id='upStatus' class='muted'></span>",
    "      </div>",
    "    </div>",
    
    # map card
    "    <div class='card'>",
    "      <h2>Location</h2>",
    paste0("      <div class='muted' style='margin-bottom:8px;'>Coordinates: ",
           ifelse(nchar(lat_txt)>0 && nchar(lon_txt)>0, paste0(lat_txt, ", ", lon_txt), "&mdash;"),
           "</div>"),
    "      <div id='map'></div>",
    "      <div id='mapStatus' class='muted' style='margin-top:8px;'></div>",
    "    </div>",
    
    "  </div>",
    
    # script
    "  <script>",
    paste0("  const API_BASE = ", jsonlite::toJSON(API_BASE, auto_unbox = TRUE), ";"),
    paste0("  const R2_PUBLIC = ", jsonlite::toJSON(R2_PUBLIC, auto_unbox = TRUE), ";"),
    paste0("  const UID = ", uid, ";"),
    paste0("  const ICE_NAME = ", jsonlite::toJSON(as.character(r$name[[1]]), auto_unbox = TRUE), ";"),
    paste0("  const ICE_LAT = ", ifelse(is.finite(lat), format(lat, scientific = FALSE), "null"), ";"),
    paste0("  const ICE_LON = ", ifelse(is.finite(lon), format(lon, scientific = FALSE), "null"), ";"),
    paste0("  const SUN_SERIES = ", sun_series_json, ";"),
    
    "  const elGal = document.getElementById('gallery');",
    "  const elImgStatus = document.getElementById('imgStatus');",
    "  const elUpStatus = document.getElementById('upStatus');",
    "  const elFile = document.getElementById('file');",
    "  const elRating = document.getElementById('rating');",
    "  const elDate = document.getElementById('shotDate');",
    "  const elComment = document.getElementById('comment');",
    "  const elFileStatus = document.getElementById('fileStatus');",
    "  const btn = document.getElementById('btnUpload');",
    "  const elSunChartHover = document.getElementById('sunChartHover');",
    
    "  function escHtml(s){",
    "    s = (s===null || s===undefined) ? '' : String(s);",
    "    return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/\"/g,'&quot;');",
    "  }",
    
    "  function ratingTxt(v){",
    "    v = Number(v);",
    "    if (v === 3) return 'good';",
    "    if (v === 2) return 'marginal';",
    "    if (v === 1) return 'poor';",
    "    return String(v);",
    "  }",
    "",
    "  function updateFileStatus(){",
    "    if (!elFileStatus) return;",
    "    const count = elFile && elFile.files ? elFile.files.length : 0;",
    "    if (!count) {",
    "      elFileStatus.textContent = 'No files selected';",
    "      return;",
    "    }",
    "    if (count === 1) {",
    "      elFileStatus.textContent = elFile.files[0].name;",
    "      return;",
    "    }",
    "    elFileStatus.textContent = `${count} files selected`;",
    "  }",
    
    "  function formatShotDate(raw){",
    "    const s = (raw === null || raw === undefined) ? '' : String(raw).trim();",
    "    if (!s) return '';",
    "    const direct = s.match(/^(\\d{4})-(\\d{2})-(\\d{2})$/);",
    "    if (direct) return `${direct[3]}.${direct[2]}.${direct[1]}`;",
    "    const parsed = new Date(s);",
    "    if (!Number.isNaN(parsed.getTime())) {",
    "      const dd = String(parsed.getDate()).padStart(2, '0');",
    "      const mm = String(parsed.getMonth() + 1).padStart(2, '0');",
    "      const yyyy = parsed.getFullYear();",
    "      return `${dd}.${mm}.${yyyy}`;",
    "    }",
    "    return s;",
    "  }",
    "",
    "  function drawSunChart(){",
    "    const canvas = document.getElementById('sunSeasonChart');",
    "    if (!canvas) return;",
    "    const ctx = canvas.getContext('2d');",
    "    if (!ctx) return;",
    "",
    "    const data = Array.isArray(SUN_SERIES) ? SUN_SERIES.slice() : [];",
    "    if (!data.length){",
    "      ctx.clearRect(0, 0, canvas.width, canvas.height);",
    "      ctx.fillStyle = '#666';",
    "      ctx.font = '14px system-ui, sans-serif';",
    "      ctx.fillText('No seasonal sun data available.', 18, 30);",
    "      if (elSunChartHover) elSunChartHover.textContent = 'No data available.';",
    "      return;",
    "    }",
    "",
    "    const w = canvas.width;",
    "    const h = canvas.height;",
    "    const pad = { t: 16, r: 20, b: 34, l: 40 };",
    "    const xMin = 0;",
    "    const xMax = data.length - 1;",
    "    const maxVal = Math.max(0.25, ...data.map(d => Number(d.sun_hours) || 0));",
    "",
    "    const xPos = (i) => {",
    "      if (xMax <= xMin) return pad.l;",
    "      return pad.l + ((i - xMin) / (xMax - xMin)) * (w - pad.l - pad.r);",
    "    };",
    "    const yPos = (v) => pad.t + (1 - (v / maxVal)) * (h - pad.t - pad.b);",
    "    const fmtMonth = (isoDate) => {",
    "      const d = new Date(String(isoDate || ''));",
    "      if (Number.isNaN(d.getTime())) return '';",
    "      return d.toLocaleDateString('de-AT', { month: 'short' });",
    "    };",
    "",
    "    const drawAxisLabels = () => {",
    "      ctx.fillStyle = '#666';",
    "      ctx.font = '11px system-ui, sans-serif';",
    "      ctx.fillText('0 h', 8, h - pad.b + 4);",
    "      ctx.fillText(`${maxVal.toFixed(1)} h`, 8, pad.t + 4);",
    "",
    "      ctx.fillStyle = '#777';",
    "      ctx.font = '10px system-ui, sans-serif';",
    "      ctx.fillText('Sun (h)', 8, h / 2);",
    "      ctx.textAlign = 'center';",
    "      ctx.fillText('Date', (pad.l + (w - pad.r)) / 2, h - 8);",
    "",
    "      const tickCount = Math.min(6, Math.max(2, data.length));",
    "      const used = new Set();",
    "      for (let i = 0; i < tickCount; i++) {",
    "        const idx = Math.round((i / (tickCount - 1 || 1)) * (data.length - 1));",
    "        if (used.has(idx)) continue;",
    "        used.add(idx);",
    "        const x = xPos(idx);",
    "        const m = fmtMonth(data[idx] && data[idx].date);",
    "        if (m) ctx.fillText(m, x, h - pad.b + 14);",
    "      }",
    "      ctx.textAlign = 'start';",
    "    };",
    "",
    "    ctx.clearRect(0, 0, w, h);",
    "",
    "    ctx.strokeStyle = '#e6e6e6';",
    "    ctx.lineWidth = 1;",
    "    for (let i = 0; i <= 4; i++){",
    "      const y = pad.t + (i / 4) * (h - pad.t - pad.b);",
    "      ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(w - pad.r, y); ctx.stroke();",
    "    }",
    "",
    "    ctx.strokeStyle = '#f1b900';",
    "    ctx.lineWidth = 2;",
    "    ctx.beginPath();",
    "    data.forEach((d, i) => {",
    "      const x = xPos(i);",
    "      const y = yPos(Number(d.sun_hours) || 0);",
    "      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);",
    "    });",
    "    ctx.stroke();",
    "",
    "    drawAxisLabels();",
    "",
    "    const drawHover = (idx) => {",
    "      ctx.clearRect(0, 0, w, h);",
    "      ctx.strokeStyle = '#e6e6e6';",
    "      ctx.lineWidth = 1;",
    "      for (let i = 0; i <= 4; i++){",
    "        const y = pad.t + (i / 4) * (h - pad.t - pad.b);",
    "        ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(w - pad.r, y); ctx.stroke();",
    "      }",
    "      ctx.strokeStyle = '#f1b900';",
    "      ctx.lineWidth = 2;",
    "      ctx.beginPath();",
    "      data.forEach((d, i) => {",
    "        const x = xPos(i);",
    "        const y = yPos(Number(d.sun_hours) || 0);",
    "        if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);",
    "      });",
    "      ctx.stroke();",
    "      drawAxisLabels();",
    "      if (idx < 0 || idx >= data.length) return;",
    "      const x = xPos(idx);",
    "      const y = yPos(Number(data[idx].sun_hours) || 0);",
    "      ctx.strokeStyle = '#999';",
    "      ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, h - pad.b); ctx.stroke();",
    "      ctx.fillStyle = '#f1b900';",
    "      ctx.beginPath(); ctx.arc(x, y, 4, 0, Math.PI * 2); ctx.fill();",
    "    };",
    "",
    "    let lastIdx = -1;",
    "    const onMove = (evt) => {",
    "      const rect = canvas.getBoundingClientRect();",
    "      const px = evt.clientX - rect.left;",
    "      const rel = Math.min(1, Math.max(0, px / rect.width));",
    "      const idx = Math.round(rel * (data.length - 1));",
    "      if (idx === lastIdx) return;",
    "      lastIdx = idx;",
    "      drawHover(idx);",
    "      const d = data[idx];",
    "      if (d && elSunChartHover) {",
    "        const hrs = Number(d.sun_hours);",
    "        elSunChartHover.textContent = `${d.date}: ${Number.isFinite(hrs) ? hrs.toFixed(2) : '0.00'} h sun`;",
    "      }",
    "    };",
    "",
    "    const onLeave = () => {",
    "      lastIdx = -1;",
    "      drawHover(-1);",
    "      if (elSunChartHover) elSunChartHover.textContent = 'Hover to show the date.';",
    "    };",
    "",
    "    drawHover(-1);",
    "    canvas.addEventListener('mousemove', onMove);",
    "    canvas.addEventListener('mouseleave', onLeave);",
    "  }",
    "",
    "  async function loadImages(){",
    "    elImgStatus.textContent = 'Loading...';",
    "    elGal.innerHTML = '';",
    "      try {",
    "      const res = await fetch(`${API_BASE}/api/images?uid=${UID}`, { method:'GET' });",
    "      if (!res.ok) throw new Error('HTTP ' + res.status);",
    "      const arr = await res.json();",
    "      if (!Array.isArray(arr) || arr.length === 0){",
    "        elImgStatus.textContent = 'No photos yet.';",
    "        return;",
    "      }",
    "      elImgStatus.textContent = `Photos: ${arr.length}`;",
    "      const groups = {};",
    "      arr.forEach(x => {",
    "        const key = [x.shot_date || '', x.rating || '', x.comment || ''].join('||');",
    "        if (!groups[key]) groups[key] = [];",
    "        groups[key].push(x);",
    "      });",
    "      const html = Object.values(groups).map(group => {",
    "        const first = group[0] || {};",
    "        const d = escHtml(formatShotDate(first.shot_date || ''));",
    "        const rt = ratingTxt(first.rating);",
    "        const c = escHtml(first.comment || '');",
    "        const cHtml = c ? `<div style=\"margin-top:6px;\">${c}</div>` : '';",
    "        const imgs = group.map(x => {",
    "          const url = x.public_url || '';",
    "          if (!url) return '';",
    "          return `",
    "            <a href=\"${url}\" target=\"_blank\" rel=\"noopener\">",
    "              <img src=\"${url}\" alt=\"Photo\" loading=\"lazy\"/>",
    "            </a>`;",
    "        }).join('');",
    "        return `",
    "          <div class=\"ph\">",
    "            <div style=\"display:grid;gap:6px;\">${imgs}</div>",
    "            <div class=\"cap\">",
    "              <b>${d}</b>",
    "              <div class=\"muted\">Condition: ${escHtml(rt)}</div>",
    "              ${cHtml}",
    "            </div>",
    "          </div>`;",
    "      }).join('');",
    "      elGal.innerHTML = html;",
    "    } catch(e){",
    "      elImgStatus.textContent = 'Loading error: ' + (e && e.message ? e.message : e);",
    "    }",
    "  }",
    
    "  function initDateDefault(){",
    "    // default = today local",
    "    const now = new Date();",
    "    const yyyy = now.getFullYear();",
    "    const mm = String(now.getMonth()+1).padStart(2,'0');",
    "    const dd = String(now.getDate()).padStart(2,'0');",
    "    if (elDate && !elDate.value) elDate.value = `${yyyy}-${mm}-${dd}`;",
    "  }",
    
    "  async function compressImage(file){",
    "    if (!file || !file.type || !file.type.startsWith('image/')) return file;",
    "    const maxDim = 1600;",
    "    const quality = 0.75;",
    "    const img = await new Promise((resolve, reject) => {",
    "      const i = new Image();",
    "      i.onload = () => resolve(i);",
    "      i.onerror = reject;",
    "      i.src = URL.createObjectURL(file);",
    "    });",
    "    const scale = Math.min(1, maxDim / Math.max(img.width, img.height));",
    "    const w = Math.max(1, Math.round(img.width * scale));",
    "    const h = Math.max(1, Math.round(img.height * scale));",
    "    const canvas = document.createElement('canvas');",
    "    canvas.width = w;",
    "    canvas.height = h;",
    "    const ctx = canvas.getContext('2d');",
    "    ctx.drawImage(img, 0, 0, w, h);",
    "    URL.revokeObjectURL(img.src);",
    "    const blob = await new Promise(resolve => canvas.toBlob(resolve, 'image/jpeg', quality));",
    "    if (!blob) return file;",
    "    const name = file.name ? file.name.replace(/\\.[^.]+$/, '') + '.jpg' : 'upload.jpg';",
    "    return new File([blob], name, { type: 'image/jpeg' });",
    "  }",
    
    "  async function doUpload(){",
    "    const files = elFile && elFile.files ? Array.from(elFile.files) : [];",
    "    if (!files.length){ elUpStatus.textContent = 'Please choose a file.'; return; }",
    "    const rating = elRating ? Number(elRating.value) : 3;",
    "    const shotDate = elDate ? String(elDate.value || '') : '';",
    "    if (!shotDate){ elUpStatus.textContent = 'Please set a date.'; return; }",
    "    const comment = (elComment && elComment.value) ? elComment.value.trim() : '';",
    
    "    let success = 0;",
    "    for (let i = 0; i < files.length; i++) {",
    "      const fRaw = files[i];",
    "    ",
    "      let f = fRaw;",
    "      try {",
    "        elUpStatus.textContent = `Image ${i+1}/${files.length}: compressing...`;",
    "        f = await compressImage(fRaw);",
    "      } catch(e) {",
    "        f = fRaw;",
    "      }",
    
    "      elUpStatus.textContent = `Image ${i+1}/${files.length}: uploading...`;",
    "      const form = new FormData();",
    "      form.append('file', f);",
    "      form.append('uid', String(UID));",
    "      form.append('ice_name', ICE_NAME || '');",
    "      form.append('rating', String(rating));",
    "      form.append('shot_date', shotDate);",
    "      if (comment) form.append('comment', comment);",
    
    "      try {",
    "        const res = await fetch(`${API_BASE}/api/upload`, { method:'POST', body: form });",
    "        const txt = await res.text();",
    "        if (!res.ok) throw new Error(txt || ('HTTP ' + res.status));",
    "        success += 1;",
    "      } catch(e){",
    "        elUpStatus.textContent = 'Upload failed: ' + (e && e.message ? e.message : e);",
    "        return;",
    "      }",
    "    }",
    "",
    "    elUpStatus.textContent = `Upload successful (${success}/${files.length}).`;",
    "    if (elFile) elFile.value = '';",
    "    if (elComment) elComment.value = '';",
    "    await loadImages();",
    "  }",
    
    "  function initMap(){",
    "    const mapStatus = document.getElementById('mapStatus');",
    "    if (ICE_LAT === null || ICE_LON === null || !isFinite(ICE_LAT) || !isFinite(ICE_LON)){",
    "      if (mapStatus) mapStatus.textContent = 'No coordinates available.';",
    "      return;",
    "    }",
    "    try {",
    "      const m = L.map('map', { scrollWheelZoom: false }).setView([ICE_LAT, ICE_LON], 13);",
    "      L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {",
    "        maxZoom: 19,",
    "        attribution: '&copy; OpenStreetMap'",
    "      }).addTo(m);",
    "      L.marker([ICE_LAT, ICE_LON]).addTo(m).bindPopup(escHtml(ICE_NAME || ('UID ' + UID))).openPopup();",
    "      if (mapStatus) mapStatus.textContent = 'Map: OSM';",
    "    } catch(e){",
    "      if (mapStatus) mapStatus.textContent = 'Map error: ' + (e && e.message ? e.message : e);",
    "    }",
    "  }",
    
    "  if (elFile) elFile.addEventListener('change', updateFileStatus);",
    "  updateFileStatus();",
    "  if (btn) btn.addEventListener('click', doUpload);",
    "  initDateDefault();",
    "  initMap();",
    "  loadImages();",
    "  drawSunChart();",
    "  </script>",
    
    "</div>",
    "</body>",
    "</html>"
  )
  
  writeLines(detail_html, out_file, useBytes = TRUE)
}

message("✅ Wrote detail pages: ", DETAIL_DIR, " (", length(unique(out$uid)), " UIDs)")



message("Done. Outputs:")
message(" - ", normalizePath(OUT_JSON, winslash = "/", mustWork = FALSE))
message(" - ", normalizePath(OUT_HTML, winslash = "/", mustWork = FALSE))
