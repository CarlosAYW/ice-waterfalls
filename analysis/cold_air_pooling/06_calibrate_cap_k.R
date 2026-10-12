# analysis/cold_air_pooling/06_calibrate_cap_k.R
# ============================================================================
# Empirical calibration of K_CAP - the Kelvin offset that corresponds to
# cap_delta_positive = 1 in cap_temperature.R.
#
# IDEA
# ----
# For any two AWS i, j the temperature difference in a well-mixed atmosphere
# is set by the lapse rate alone:
#
#     T_i - T_j  =  -gamma * (z_i - z_j)
#
# The residual of that relation is what the lapse rate cannot explain:
#
#     r_ij  =  (T_i - T_j) + gamma * (z_i - z_j)
#
# The station sitting deeper in the cold pool is colder than the lapse rate
# predicts, so r should scale with the difference in topographic CAP
# potential - but only to the degree that a CAP is actually developed at that
# hour. The model in cap_temperature.R writes that as
#
#     dT_cap = -K_CAP * cap_delta * gate,     gate = f_inv * f_wind * f_ramp
#
# so the estimator has to use the SAME gate, otherwise K_CAP is attenuated by
# the mean gate value and then attenuated a second time when the model applies
# the gate again. A synthetic test with K_CAP = 6.8 K recovers 1.80 K from the
# naive regression on dcap alone, and 6.78 K from the gated regression:
#
#     r_ij  =  -K_CAP * (dcap * gate)  +  eps          <- regress on THIS
#
# Both fits are reported so the attenuation is visible and quotable.
#
# Two controls make the result trustworthy:
#   (a) gamma is FITTED from well-mixed (windy, non-inversion) periods rather
#       than assumed, and
#   (b) the same regression is run on well-mixed periods, where the slope must
#       be indistinguishable from zero.
#
# OUTPUT
# ------
#   data/CAP/cap_k_calibration.csv        one row per subset with slope + CI
#   data/CAP/cap_k_calibration_pairs.csv  the pair-level sample actually used
#   output/cap/cap_k_calibration.png      diagnostic scatter + fit
#
# Run from the repository root:
#   Rscript analysis/cold_air_pooling/06_calibrate_cap_k.R
#   Rscript analysis/cold_air_pooling/06_calibrate_cap_k.R --seasons=5 --force
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(purrr)
  library(httr2)
})

source("analysis/cold_air_pooling/cap_temperature.R")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

TZ_LOCAL <- "Europe/Vienna"

PATH_STATIONS  <- "data/AWS/stations_all.csv"
PATH_CAP_INDEX <- "data/CAP/cap_index.csv"
PATH_ICE_PAIRS <- "data/CAP/icefall_station_cap.csv"
DIR_INV_CACHE  <- "data/_cache_inversion"
DIR_AWS_CACHE  <- "data/_cache_aws_tl"      # created by this script
OUT_DIR        <- "data/CAP"
OUT_PLOT_DIR   <- "output/cap"

MAX_PAIR_DIST_KM <- 30    # keep pairs under comparable synoptic conditions
MIN_PAIR_DZ_M    <- 100   # need height contrast to constrain gamma
MAX_PAIR_DZ_M    <- 1600  # avoid free-atmosphere / summit stations
MIN_ABS_DCAP     <- 0.05  # near-identical CAP scores carry no signal
MAX_STATIONS     <- 40    # cap the API load

WINTER_MONTHS <- c(11, 12, 1, 2, 3)
NIGHT_HOURS   <- c(19:23, 0:6)      # local time
N_SEASONS     <- 3

CALM_FF_MS       <- 1.0   # Rauchoecker et al. 2024: CAP undisturbed
WINDY_FF_MS      <- 3.0   # well-mixed control subset
MIN_OBS_PER_FIT  <- 200
N_BOOT           <- 500

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
ns_arg <- grep("^--seasons=", args, value = TRUE)
if (length(ns_arg)) N_SEASONS <- as.integer(sub("^--seasons=", "", ns_arg[[1]]))
force_refetch <- any(args == "--force")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_AWS_CACHE, recursive = TRUE, showWarnings = FALSE)

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371; p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 +
    cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

# ---------------------------------------------------------------------------
# 1) Station selection
# ---------------------------------------------------------------------------

stations_all <- read_csv(PATH_STATIONS, show_col_types = FALSE, progress = FALSE)
cap_index    <- read_csv(PATH_CAP_INDEX, show_col_types = FALSE, progress = FALSE)

station_cap <- cap_index %>%
  filter(entity_type == "station", qa_flag == "ok") %>%
  transmute(station_id = as.character(entity_id), cap_potential, cap_confidence)

# Only GeoSphere stations are used: the LWD Tirol reader lives in
# scripts/00_build_inversion_cache.R and is not available here, so LWD stations
# would silently yield no data. This costs coverage (49 of the 75 stations with
# a CAP score are LWD) and must be stated as a limitation of the calibration.
stations <- stations_all %>%
  mutate(station_id = as.character(station_id),
         has_tl = as.logical(tl),
         has_ff = as.logical(ff)) %>%
  filter(!is.na(has_tl), has_tl) %>%
  inner_join(station_cap, by = "station_id") %>%
  filter(is.finite(altitude_m), is.finite(lat), is.finite(lon),
         is.finite(cap_potential))

n_lwd <- sum(stations$source != "GeoSphere", na.rm = TRUE)
stations <- stations %>%
  filter(source == "GeoSphere") %>%
  arrange(desc(coalesce(has_ff, FALSE)), desc(cap_confidence)) %>%
  slice_head(n = MAX_STATIONS)

message("Calibration stations: ", nrow(stations),
        " GeoSphere (", sum(stations$has_ff, na.rm = TRUE), " with wind); ",
        n_lwd, " LWD stations skipped (no reader available)")
if (nrow(stations) < 5) {
  stop("Need >= 5 GeoSphere stations with a CAP score. Run 05_build_cap_index.R ",
       "with --all-stations first.")
}

# ---------------------------------------------------------------------------
# 2) Fetch / cache hourly TL (+ FF) per station
# ---------------------------------------------------------------------------

END_DATE <- as.Date(with_tz(Sys.time(), TZ_LOCAL))
start_year <- year(END_DATE) - N_SEASONS - as.integer(month(END_DATE) < 10)
START_DATE <- as.Date(sprintf("%d-10-01", start_year))
message("Period: ", START_DATE, " .. ", END_DATE)

fetch_station_series <- function(station_id, source = "GeoSphere") {
  cache_file <- file.path(
    DIR_AWS_CACHE,
    sprintf("%s_%s_%s.rds", gsub("[^A-Za-z0-9]", "", station_id),
            format(START_DATE, "%Y%m%d"), format(END_DATE, "%Y%m%d")))
  if (file.exists(cache_file) && !force_refetch) return(readRDS(cache_file))

  # Same dataset and query shape as scripts/00_build_inversion_cache.R, which is
  # known to work against the GeoSphere API in this project.
  # Ask for tl+ff; if the station does not serve ff the whole request fails, so
  # fall back to tl only (the wind gate then uses its NA default).
  geosphere_csv <- function(pars) {
    req <- request("https://dataset.api.hub.geosphere.at/v1/station/historical/klima-v2-10min") |>
      req_url_query(station_ids = as.character(station_id),
                    parameters = pars,
                    start = as.character(START_DATE),
                    end = as.character(END_DATE + 1),
                    output_format = "csv") |>
      req_user_agent("icefall-model/1.0 (R httr2)") |>
      req_retry(max_tries = 3)
    resp <- req_perform(req)
    if (resp_status(resp) >= 400) {
      stop("HTTP ", resp_status(resp), ": ", substr(resp_body_string(resp), 1, 200))
    }
    readr::read_csv(I(resp_body_string(resp)), show_col_types = FALSE,
                    progress = FALSE)
  }

  out <- tryCatch({
    df <- tryCatch(geosphere_csv("tl,ff"), error = function(e) {
      message("  ~ ", station_id, ": tl,ff failed (", conditionMessage(e),
              ") - retrying with tl only")
      geosphere_csv("tl")
    })
    nms <- tolower(names(df))
    tcol <- names(df)[which.max(grepl("time|timestamp|date", nms))]
    i_tl <- which(nms == "tl")[1]; i_ff <- which(nms == "ff")[1]
    if (is.na(i_tl)) stop("no 'tl' column; got: ", paste(names(df), collapse = ", "))
    tibble(
      timestamp = as.POSIXct(df[[tcol]], tz = "UTC"),
      TL = suppressWarnings(as.numeric(df[[i_tl]])),
      FF = if (!is.na(i_ff)) suppressWarnings(as.numeric(df[[i_ff]])) else NA_real_)
  }, error = function(e) {
    message("  ! ", station_id, " FAILED: ", conditionMessage(e)); NULL
  })

  if (is.null(out) || nrow(out) == 0) {
    message("  ! ", station_id, ": 0 rows returned")
    return(NULL)
  }
  out <- out %>%
    mutate(timestamp = with_tz(timestamp, TZ_LOCAL),
           time = floor_date(timestamp, "hour")) %>%
    group_by(time) %>%
    summarise(TL = mean(TL, na.rm = TRUE), FF = mean(FF, na.rm = TRUE),
              .groups = "drop") %>%
    filter(is.finite(TL))
  saveRDS(out, cache_file)
  out
}

message("Fetching station series (", nrow(stations), " stations) ...")
n_ok <- 0
series <- map2(stations$station_id, stations$source, function(sid, src) {
  s <- fetch_station_series(sid, src)
  if (is.null(s) || nrow(s) == 0) return(NULL)
  n_ok <<- n_ok + 1
  message(sprintf("  ok %-10s %6d hours  [%s .. %s]", sid, nrow(s),
                  format(min(s$time), "%Y-%m-%d"), format(max(s$time), "%Y-%m-%d")))
  s$station_id <- sid
  s
}) %>% compact() %>% bind_rows()

if (nrow(series) == 0) {
  stop("No AWS data retrieved from any of the ", nrow(stations), " stations.\n",
       "Check the messages above: an HTTP error points at the API/endpoint, ",
       "'0 rows returned' at the requested period. A quick manual probe:\n",
       "  httr2::request(\"https://dataset.api.hub.geosphere.at/v1/station/historical/klima-v2-10min\") |>\n",
       "    httr2::req_url_query(station_ids=\"38\", parameters=\"tl\", start=\"2026-01-01\", ",
       "end=\"2026-01-02\", output_format=\"csv\") |> httr2::req_perform() |> httr2::resp_body_string()")
}
message("  rows: ", nrow(series), " for ", n_distinct(series$station_id), " stations")
stations <- stations %>% filter(station_id %in% unique(series$station_id))

# ---------------------------------------------------------------------------
# 3) Inversion state + hours since CAP onset
# ---------------------------------------------------------------------------

inv_files <- list.files(DIR_INV_CACHE, pattern = "^inversion_\\d{8}\\.rds$",
                        full.names = TRUE)
if (length(inv_files) == 0) {
  stop("No inversion cache in ", DIR_INV_CACHE,
       ". Run scripts/00_build_inversion_cache.R first.")
}
inv <- readRDS(inv_files[which.max(file.mtime(inv_files))]) %>%
  mutate(time = floor_date(time, "hour")) %>%
  group_by(time) %>%
  summarise(inv_active = any(inv_active, na.rm = TRUE),
            inv_grad   = mean(inv_grad_max_C_per_100m, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(time)

# hours since inv_active turned on / off (hourly grid -> step_min = 60)
clock <- cap_state_clock(inv$time, inv$inv_active, step_min = 60)
inv <- bind_cols(inv, clock[, c("since_onset_h", "since_erosion_h")])
message("Inversion cache: ", format(min(inv$time)), " .. ", format(max(inv$time)))

# ---------------------------------------------------------------------------
# 4) Pair sample + the model gate
# ---------------------------------------------------------------------------

# The gate depends on inv_active/inv_grad, which only exist for the period the
# inversion cache covers. Rows outside it would silently get gate = 0 and only
# dilute the regression, so they are dropped here.
inv_from <- min(inv$time); inv_to <- max(inv$time)
n_before <- nrow(series)
series <- series %>% filter(time >= inv_from, time <= inv_to)
message(sprintf("Restricted to inversion-cache period %s .. %s: %d of %d rows kept",
                format(inv_from, "%Y-%m-%d"), format(inv_to, "%Y-%m-%d"),
                nrow(series), n_before))
if (nrow(series) == 0) {
  stop("No station data inside the inversion-cache period. Rebuild the cache ",
       "for the seasons you want to calibrate on, or lower --seasons.")
}

obs <- series %>%
  left_join(inv, by = "time") %>%
  mutate(month = month(time), hour = hour(time)) %>%
  filter(month %in% WINTER_MONTHS, hour %in% NIGHT_HOURS)

st_small <- stations %>%
  transmute(station_id, lat, lon, z = altitude_m, cap = cap_potential)

pairs_meta <- crossing(i = st_small$station_id, j = st_small$station_id) %>%
  filter(i < j) %>%
  left_join(st_small %>% rename_with(~paste0(.x, "_i")), by = c("i" = "station_id_i")) %>%
  left_join(st_small %>% rename_with(~paste0(.x, "_j")), by = c("j" = "station_id_j")) %>%
  mutate(dist_km = haversine_km(lat_i, lon_i, lat_j, lon_j),
         dz_m = z_i - z_j,
         dcap = cap_i - cap_j) %>%
  filter(dist_km <= MAX_PAIR_DIST_KM,
         abs(dz_m) >= MIN_PAIR_DZ_M, abs(dz_m) <= MAX_PAIR_DZ_M,
         abs(dcap) >= MIN_ABS_DCAP) %>%
  transmute(station_id_i = i, station_id_j = j, dist_km, dz_m, dcap)

message("Candidate station pairs: ", nrow(pairs_meta))
if (nrow(pairs_meta) == 0) {
  stop("No usable station pairs. Relax MAX_PAIR_DIST_KM / MIN_ABS_DCAP.")
}

sample_df <- pairs_meta %>%
  inner_join(obs %>% select(time, station_id, TL, FF, inv_active, inv_grad,
                            since_onset_h, since_erosion_h),
             by = c("station_id_i" = "station_id"),
             relationship = "many-to-many") %>%
  rename(TL_i = TL, FF_i = FF) %>%
  inner_join(obs %>% select(time, station_id, TL, FF),
             by = c("time", "station_id_j" = "station_id")) %>%
  rename(TL_j = TL, FF_j = FF) %>%
  mutate(
    dT_obs     = TL_i - TL_j,
    ff_max     = pmax(FF_i, FF_j, na.rm = TRUE),
    ff_max     = ifelse(is.finite(ff_max), ff_max, NA_real_),
    inv_active = coalesce(inv_active, FALSE),
    inv_grad   = coalesce(inv_grad, 0)
  ) %>%
  filter(is.finite(dT_obs), is.finite(dz_m), is.finite(dcap))

# the SAME gate the model applies - this is what makes the fit self-consistent
sample_df <- sample_df %>%
  mutate(
    f_inv  = cap_inversion_factor(inv_grad, inv_active, mode = "hybrid"),
    f_wind = cap_wind_factor(ff_max),
    f_ramp = cap_ramp_factor(since_onset_h, since_erosion_h),
    gate   = f_inv * f_wind * f_ramp,
    dcap_gated = dcap * gate,
    # second regressor: lets the effective lapse rate respond to stability
    dz_gate    = dz_m * gate
  )

message("Pair-hour observations: ", nrow(sample_df),
        " | mean gate = ", sprintf("%.3f", mean(sample_df$gate, na.rm = TRUE)))

# ---------------------------------------------------------------------------
# 5) gamma from WELL-MIXED conditions (fitted, not assumed)
# ---------------------------------------------------------------------------

wellmixed <- sample_df %>%
  filter(!inv_active, inv_grad < CAP_PARAMS$INV_GRAD_OFF,
         is.finite(ff_max), ff_max >= WINDY_FF_MS)

if (nrow(wellmixed) >= MIN_OBS_PER_FIT) {
  GAMMA_FIT <- -unname(coef(lm(dT_obs ~ 0 + dz_m, data = wellmixed))[1])
  message(sprintf("gamma fitted from %d well-mixed obs: %.5f K/m (%.2f K/100 m)",
                  nrow(wellmixed), GAMMA_FIT, GAMMA_FIT * 100))
} else {
  GAMMA_FIT <- 0.0065
  message("Too few well-mixed obs (", nrow(wellmixed),
          ") - falling back to gamma = 0.0065 K/m")
}
if (!is.finite(GAMMA_FIT) || GAMMA_FIT < 0.002 || GAMMA_FIT > 0.012) {
  warning(sprintf("Fitted gamma = %.5f K/m is implausible; using 0.0065.", GAMMA_FIT))
  GAMMA_FIT <- 0.0065
}

# Two residuals. resid_K uses the lapse rate fitted from windy, inversion-free
# periods; resid_K_std uses the standard 0.0065 K/m. In Alpine valleys even
# "well-mixed" winter nights retain some stability, so the fitted gamma tends to
# be too shallow and absorbs part of the cold-air signal - which biases K_CAP
# low. Reporting both makes that bias visible instead of hiding it.
GAMMA_STD <- 0.0065
sample_df <- sample_df %>%
  mutate(resid_K     = dT_obs + GAMMA_FIT * dz_m,
         resid_K_std = dT_obs + GAMMA_STD * dz_m)

# rebuild wellmixed from the updated frame so it carries the residual columns
wellmixed <- sample_df %>%
  filter(!inv_active, inv_grad < CAP_PARAMS$INV_GRAD_OFF,
         is.finite(ff_max), ff_max >= WINDY_FF_MS)

# ---------------------------------------------------------------------------
# 6) Regressions
# ---------------------------------------------------------------------------

# Blocked bootstrap over station PAIRS, not over individual pair-hours.
# The 53 pairs observe the same 1811 night hours, so pair-hours are strongly
# dependent. Resampling rows treats them as independent and produces intervals
# that are far too narrow (5.87-6.76 K instead of 1.6-14.4 K for the same fit).
boot_ci <- function(d, xcols, ycol = "resid_K", n = N_BOOT) {
  if (nrow(d) < MIN_OBS_PER_FIT) return(c(NA_real_, NA_real_))
  f <- as.formula(paste(ycol, "~ 0 +", paste(xcols, collapse = " + ")))
  blocks <- split(seq_len(nrow(d)), paste(d$station_id_i, d$station_id_j))
  nb <- length(blocks)
  s <- replicate(n, {
    idx <- unlist(blocks[sample.int(nb, nb, replace = TRUE)], use.names = FALSE)
    cf <- try(coef(lm(f, data = d[idx, , drop = FALSE])), silent = TRUE)
    if (inherits(cf, "try-error")) NA_real_ else unname(cf[1])
  })
  unname(quantile(s, c(0.025, 0.975), na.rm = TRUE))
}

# A single fixed lapse rate cannot describe both well-mixed and strongly stable
# nights. Whatever it gets wrong scales with dz and lands in the residual, and
# because dz and dcap are correlated across pairs (r = -0.27) it leaks into the
# CAP slope. Adding dz_gate lets the effective gradient steepen with stability
# and removes that leak: a synthetic test with a known K_CAP = 6.5 K recovers
# 2.7 K from the one-term fit on calm nights but 6.6 K from the two-term fit.
fit_subset <- function(d, label, xcol = "dcap_gated", ycol = "resid_K",
                       with_dz = FALSE) {
  xcols <- if (with_dz) c(xcol, "dz_gate") else xcol
  d <- d %>% filter(if_all(all_of(c(ycol, xcols)), is.finite))
  gam <- if (ycol == "resid_K") GAMMA_FIT else GAMMA_STD
  base <- tibble(subset = label, regressor = xcol, lapse_ref = ycol,
                 model = if (with_dz) "gated + dz_gate" else "gated only",
                 gamma_used_K_per_m = gam, n_obs = nrow(d),
                 n_pairs = n_distinct(paste(d$station_id_i, d$station_id_j)))
  if (nrow(d) < MIN_OBS_PER_FIT ||
      any(vapply(xcols, function(c_) sd(d[[c_]], na.rm = TRUE) == 0, logical(1)))) {
    return(bind_cols(base, tibble(K_CAP = NA_real_, ci_lo = NA_real_,
                                  ci_hi = NA_real_, dgamma_K_per_100m = NA_real_,
                                  r2 = NA_real_,
                                  mean_gate = mean(d$gate, na.rm = TRUE))))
  }
  f  <- lm(as.formula(paste(ycol, "~ 0 +", paste(xcols, collapse = " + "))), data = d)
  ci <- boot_ci(d, xcols, ycol)
  bind_cols(base, tibble(
    K_CAP = -unname(coef(f)[1]),
    ci_lo = -ci[2], ci_hi = -ci[1],
    dgamma_K_per_100m = if (with_dz) unname(coef(f)[2]) * 100 else NA_real_,
    r2 = summary(f)$r.squared,
    mean_gate = mean(d$gate, na.rm = TRUE)))
}

subsets <- list(
  "ALL winter nights (gated)"              = sample_df,
  "CAP nights (inv_active)"                = sample_df %>% filter(inv_active),
  "CAP nights, calm (<=1 m/s)"             = sample_df %>% filter(inv_active, ff_max <= CALM_FF_MS),
  "strong inversion (>=0.5 C/100m), calm"  = sample_df %>% filter(inv_grad >= 0.50, ff_max <= CALM_FF_MS),
  "moderate inversion (0.2-0.5), calm"     = sample_df %>% filter(inv_grad >= 0.20, inv_grad < 0.50, ff_max <= CALM_FF_MS)
)

# c() first: bind_rows() must receive ONE flat list of tibbles. Passing several
# lists as separate arguments makes it treat each list as a single row and
# produces list-columns (which then break write_csv).
fit_list <- c(
  imap(subsets, ~fit_subset(.x, .y, "dcap_gated", "resid_K")),
  # the same fits with the stability-dependent lapse-rate term - these are the
  # ones reported in the thesis, because they are subset-independent
  imap(subsets, ~fit_subset(.x, .y, "dcap_gated", "resid_K", with_dz = TRUE)),
  imap(subsets, ~fit_subset(.x, .y, "dcap_gated", "resid_K_std")),
  # naive (ungated) fit for comparison - shows the attenuation
  imap(subsets, ~fit_subset(.x, .y, "dcap", "resid_K")),
  # control: well-mixed. gate is 0 here, so both regressors would be identically
  # zero; the ungated quantities are used instead. Must still return ~0.
  list(fit_subset(wellmixed %>% mutate(dcap_gated = dcap, dz_gate = dz_m),
                  "CONTROL: well-mixed & windy", "dcap_gated", "resid_K")),
  list(fit_subset(wellmixed %>% mutate(dcap_gated = dcap, dz_gate = dz_m),
                  "CONTROL: well-mixed & windy", "dcap_gated", "resid_K",
                  with_dz = TRUE))
)
stopifnot(all(vapply(fit_list, is.data.frame, logical(1))))
results <- bind_rows(fit_list)
stopifnot(!any(vapply(results, is.list, logical(1))))

results <- results %>%
  mutate(
    gamma_fit_K_per_m = GAMMA_FIT,
    lit_lo = CAP_PARAMS$K_CAP_MIN_LIT,
    lit_hi = CAP_PARAMS$K_CAP_MAX_LIT,
    within_literature = is.finite(K_CAP) &
      K_CAP >= CAP_PARAMS$K_CAP_MIN_LIT & K_CAP <= CAP_PARAMS$K_CAP_MAX_LIT,
    cap_algorithm_version = "cap_dem_v1"
  )

show <- results %>%
  select(subset, model, regressor, lapse_ref, n_obs, n_pairs, mean_gate,
         K_CAP, ci_lo, ci_hi, dgamma_K_per_100m, r2)
print(as.data.frame(show), row.names = FALSE, digits = 3)

n_skipped <- sum(is.na(results$K_CAP))
if (n_skipped > 0) {
  message("\n", n_skipped, " of ", nrow(results), " fits were skipped for having ",
          "fewer than ", MIN_OBS_PER_FIT, " observations - see n_obs above.")
}

if (GAMMA_FIT < 0.0045) {
  message("\nNOTE: the fitted lapse rate (", sprintf("%.2f", GAMMA_FIT * 100),
          " K/100 m) is well below the free-atmosphere value of 0.65 K/100 m. ",
          "This is typical for Alpine valleys in winter, but it means the ",
          "'well-mixed' reference still contains cold-air pooling, so the K_CAP ",
          "fitted against it is a LOWER bound. Compare the lapse_ref = ",
          "'resid_K_std' rows, which use 0.0065 K/m instead.")
}

write_csv(results, file.path(OUT_DIR, "cap_k_calibration.csv"))
write_csv(sample_df %>% slice_sample(n = min(nrow(sample_df), 200000)),
          file.path(OUT_DIR, "cap_k_calibration_pairs.csv"))

# ---------------------------------------------------------------------------
# 7) Diagnostic plot
# ---------------------------------------------------------------------------

MAIN <- "CAP nights, calm (<=1 m/s)"
plot_dat <- subsets[[MAIN]] %>% filter(is.finite(resid_K), is.finite(dcap_gated))
sel <- results$subset == MAIN & results$regressor == "dcap_gated" &
       results$lapse_ref == "resid_K"
K_fit <- results$K_CAP[sel][1]

if (nrow(plot_dat) >= MIN_OBS_PER_FIT && is.finite(K_fit)) {
  p <- ggplot(plot_dat, aes(dcap_gated, resid_K)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    annotate("ribbon", x = c(-1, 1),
             ymin = c(CAP_PARAMS$K_CAP_MAX_LIT, -CAP_PARAMS$K_CAP_MAX_LIT),
             ymax = c(CAP_PARAMS$K_CAP_MIN_LIT, -CAP_PARAMS$K_CAP_MIN_LIT),
             alpha = 0.10, fill = "orange") +
    geom_point(alpha = 0.06, size = 0.6, colour = "steelblue4") +
    geom_smooth(method = "lm", formula = y ~ 0 + x, se = TRUE,
                colour = "firebrick", linewidth = 0.9) +
    labs(
      title = "Empirical calibration of the CAP temperature scaling K_CAP",
      subtitle = sprintf("Calm inversion nights (Nov-Mar, night hours). K_CAP = %.2f K [%.2f, %.2f]; literature 5-9 K",
                         K_fit,
                         results$ci_lo[sel][1], results$ci_hi[sel][1]),
      x = expression(Delta*"CAP potential "%*%" gate  ("*f[inv]%.%f[wind]%.%f[ramp]*")"),
      y = expression("Lapse-rate residual  "*r[ij]*"  (K)"),
      caption = sprintf(paste0("gamma = %.2f K/100 m fitted from well-mixed periods; ",
                               "n = %d pair-hours; orange band = literature range ",
                               "(Rauchoecker et al. 2024)"),
                        GAMMA_FIT * 100, nrow(plot_dat))
    ) +
    theme_minimal(base_size = 11)
  ggsave(file.path(OUT_PLOT_DIR, "cap_k_calibration.png"), p,
         width = 8, height = 5.2, dpi = 200)
  message("Wrote ", file.path(OUT_PLOT_DIR, "cap_k_calibration.png"))
}

# ---------------------------------------------------------------------------
# 8) What the fitted K_CAP means for the actual icefalls
# ---------------------------------------------------------------------------

if (is.finite(K_fit) && file.exists(PATH_ICE_PAIRS)) {
  d <- read_csv(PATH_ICE_PAIRS, show_col_types = FALSE, progress = FALSE)$cap_delta_positive
  d <- d[is.finite(d)]
  message("\n--- Implied CAP correction per icefall (fully developed CAP, gate = 1) ---")
  message(sprintf("  K_CAP (fitted) = %.2f K   [literature 5.0-9.0 K]", K_fit))
  message(sprintf("  median %.2f K | p75 %.2f K | p90 %.2f K | max %.2f K",
                  K_fit * median(d), K_fit * quantile(d, .75),
                  K_fit * quantile(d, .90), K_fit * max(d)))
  message(sprintf("  icefalls with |dT| > 1 K: %d / %d", sum(d * K_fit > 1), length(d)))

  naive_K <- results$K_CAP[results$subset == MAIN & results$regressor == "dcap" &
                             results$lapse_ref == "resid_K"][1]
  if (is.finite(naive_K)) {
    message(sprintf("\n  (ungated fit would give %.2f K - the attenuation by the mean\n",
                    naive_K),
            sprintf("   gate of %.3f is why the gated regressor must be used)",
                    mean(plot_dat$gate, na.rm = TRUE)))
  }

  if (!isTRUE(results$within_literature[sel][1])) {
    warning("Fitted K_CAP is outside the literature band - discuss why ",
            "(sparse station network, few CAP nights, DEM index bias, ",
            "stations not representative of icefall micro-sites).")
  }
  ctrl <- results$K_CAP[results$subset == "CONTROL: well-mixed & windy"][1]
  if (is.finite(ctrl) && abs(ctrl) > 0.25 * abs(K_fit)) {
    warning(sprintf(paste0("Control slope is %.2f K, i.e. >25%% of the CAP-night ",
                           "slope. Part of the signal is NOT cold-air pooling ",
                           "(residual lapse-rate or station-exposure bias)."), ctrl))
  }
}

gated <- results %>%
  filter(regressor == "dcap_gated", subset != "CONTROL: well-mixed & windy",
         is.finite(K_CAP))
if (nrow(gated) > 1) {
  message("\n--- Spread across defensible specifications ---")
  message(sprintf("  K_CAP ranges from %.2f to %.2f K over %d subset/lapse-rate combinations.",
                  min(gated$K_CAP), max(gated$K_CAP), nrow(gated)))
  message("  A range this wide means the fit is NOT a point estimate. Report the ",
          "range and run the model at both ends rather than quoting one value.")
}

message("\nWrote ", file.path(OUT_DIR, "cap_k_calibration.csv"))
message("Set CAP_PARAMS$K_CAP in cap_temperature.R to the fitted value, or pass ",
        "K_CAP = <value> to cap_temperature_offset().")
