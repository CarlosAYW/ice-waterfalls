# analysis/cold_air_pooling/cap_temperature.R
# ============================================================================
# Convert the dimensionless, purely topographic CAP potential (cap_dem_v1,
# see 05_build_cap_index.R) into a PHYSICAL temperature offset in Kelvin.
#
# Until now the icefall temperature was derived from its assigned AWS by a
# constant lapse rate only:
#
#     T_icefall = T_station - GAMMA * (z_icefall - z_station)
#
# This module adds the cold-air-pool (CAP) term:
#
#     T_icefall = T_station - GAMMA * dz + dT_cap
#     dT_cap    = -K_CAP * cap_delta_positive * f_inv * f_wind * f_ramp
#                 * f_height * f_conf                                   [K]
#
# All coefficients are anchored in measured Alpine CAP literature so that the
# scaling is defensible rather than tuned:
#
#   Rauchoecker, A., Rudolph, A., Stiperski, I. & Lehner, M. (2024):
#     Cold-air pool development in a small Alpine valley.
#     Q. J. R. Meteorol. Soc. 150, 1243-1266. doi:10.1002/qj.4644
#     -> Seefeld (Tyrol), valley 200-300 m deep, floor ~1180-1270 m a.s.l.
#        * dT(valley floor vs. slope stations) = 5-9 K
#        * initial cooling rate 3.0 K/h, cooling time constant tau = 2.9 h
#          (Gruenloch 3.2 h, Peter Sinks 3.1 h)
#        * near-surface gradient inside the inversion ~2 degC/10 m
#        * inversion depth ~30 m, isothermal above
#        * wind > 1 m/s already warms the pool, > 2 m/s erodes it entirely
#
#   Pospichal, B., Eisenbach, S., Whiteman, C.D., Steinacker, R. &
#     Dorninger, M. (2003): Observations of the cold air outflow from a basin
#     cold pool through a low pass. (Gruenloch / Gstettneralm sinkhole, Lunz)
#     -> potential-temperature inversions > 25 K over 80-120 m depth.
#        Used here ONLY as the upper physical bound for extreme, doline-like,
#        strongly sheltered sites - NOT as the default scaling.
#
# Author: Carlos Wydra - MSc thesis, Univ. Graz
# ============================================================================

# ---------------------------------------------------------------------------
# 1) Literature-anchored constants
# ---------------------------------------------------------------------------

CAP_PARAMS <- list(

  # Temperature offset [K] that corresponds to cap_delta_positive = 1.
  # Seefeld measured 5-9 K between valley floor and slope in a 200-300 m deep
  # Alpine valley -> midpoint 6.5 K is the default for ordinary valley sites.
  # Overwrite with the empirically fitted value from 06_calibrate_cap_k.R.
  K_CAP = 6.5,

  # Plausibility band for K_CAP. The upper bound is the Gruenloch extreme and
  # should only ever be approached by doline-like, fully enclosed basins.
  K_CAP_MIN_LIT = 5.0,    # Rauchoecker et al. 2024, lower end
  K_CAP_MAX_LIT = 9.0,    # Rauchoecker et al. 2024, upper end
  K_CAP_HARD_MAX = 25.0,  # Pospichal et al. 2003 (Gruenloch), absolute ceiling

  # --- Inversion gating (uses the global inversion cache) -------------------
  # Thresholds are the class boundaries already used in
  # 00_build_inversion_cache.R (bulk gradient Imst 860 m -> Muttekopfhuette
  # 1935 m, i.e. over dz = 1075 m). NOTE: this is a BULK gradient over the
  # whole valley relief and must not be compared directly with the
  # near-surface values in the literature (see the thesis Methods section).
  INV_GRAD_OFF = 0.05,    # degC/100 m -> "none"
  INV_GRAD_FULL = 0.50,   # degC/100 m -> "strong"

  # --- Wind gating (Rauchoecker et al. 2024) --------------------------------
  WIND_FULL_MS = 1.0,     # <= 1 m/s: CAP undisturbed
  WIND_ZERO_MS = 2.0,     # >= 2 m/s: CAP eroded

  # --- Temporal build-up ----------------------------------------------------
  TAU_BUILD_H = 2.9,      # cooling time constant, Seefeld
  TAU_ERODE_H = 1.0,      # faster erosion after sunrise / wind onset

  # --- Vertical structure (OPTIONAL refinement, off by default) -------------
  # CAUTION - scale mismatch. CAP_DEPTH_M = 30 m is the NEAR-SURFACE inversion
  # depth measured at Seefeld (~2 degC/10 m inside, isothermal above), whereas
  # cap_delta_positive is derived from a DEM window of 1400 m radius and
  # therefore describes pooling on the VALLEY cross-section scale. Multiplying
  # the two would discount the same signal twice. Height weighting is thus
  # switched OFF in the default correction and is provided as a separate
  # sensitivity experiment for the vertical-structure discussion.
  CAP_DEPTH_M = 30,       # Rauchoecker et al. 2024, near-surface inversion depth
  DEFAULT_ICEFALL_HEIGHT_M = 30,
  DEFAULT_ICEFALL_BASE_M = 0    # base height above the local valley floor
)

clamp01 <- function(x) pmax(0, pmin(1, x))


# ---------------------------------------------------------------------------
# 2) Individual gating factors (all return values in [0, 1])
# ---------------------------------------------------------------------------

#' Inversion factor from the global inversion cache.
#'
#' @param inv_grad_C_per_100m bulk inversion gradient (degC/100 m)
#' @param inv_active logical, hysteresis-filtered inversion state
#' @param mode "hybrid" (default): continuous ramp, but zero unless inv_active;
#'   "continuous": ramp only; "binary": inv_active only.
cap_inversion_factor <- function(inv_grad_C_per_100m,
                                 inv_active = NULL,
                                 mode = c("hybrid", "continuous", "binary"),
                                 p = CAP_PARAMS) {
  mode <- match.arg(mode)
  g <- as.numeric(inv_grad_C_per_100m)
  ramp <- clamp01((g - p$INV_GRAD_OFF) / (p$INV_GRAD_FULL - p$INV_GRAD_OFF))
  ramp[!is.finite(ramp)] <- 0

  if (mode == "continuous" || is.null(inv_active)) return(ramp)

  act <- as.logical(inv_active)
  act[is.na(act)] <- FALSE
  if (mode == "binary") return(as.numeric(act))
  ramp * as.numeric(act)
}


#' Wind factor. 1 below WIND_FULL_MS, linearly to 0 at WIND_ZERO_MS.
#'
#' Rauchoecker et al. (2024) observed warming in the pool already above
#' ~1 m/s and complete erosion of the CAP above ~2 m/s (foehn breakthrough).
#'
#' @param ff_ms wind speed [m/s]. NA -> na_value (default 0.5, i.e. a
#'   conservative half-strength CAP when no wind measurement exists).
cap_wind_factor <- function(ff_ms, na_value = 0.5, p = CAP_PARAMS) {
  f <- clamp01((p$WIND_ZERO_MS - as.numeric(ff_ms)) /
                 (p$WIND_ZERO_MS - p$WIND_FULL_MS))
  f[!is.finite(f)] <- na_value
  f
}


#' Exponential build-up / erosion of the pool.
#'
#' Default driver is the time since inv_active switched on (needs no
#' astronomy). Set `since_sunset_h` instead to drive it from sunset.
#'
#' @param since_onset_h hours since the CAP-favourable state began
#' @param since_erosion_h hours since the state ended (NA while still active)
cap_ramp_factor <- function(since_onset_h,
                            since_erosion_h = NA_real_,
                            p = CAP_PARAMS) {
  build <- 1 - exp(-pmax(0, as.numeric(since_onset_h)) / p$TAU_BUILD_H)
  build[!is.finite(build)] <- 0

  erode <- exp(-pmax(0, as.numeric(since_erosion_h)) / p$TAU_ERODE_H)
  erode[!is.finite(erode)] <- 1

  clamp01(build * erode)
}


#' Helper: hours since a logical state turned TRUE / FALSE, for a regular
#' time series. Returns a data.frame with since_onset_h / since_erosion_h.
cap_state_clock <- function(time, state, step_min = 10) {
  state <- as.logical(state); state[is.na(state)] <- FALSE
  n <- length(state)
  since_on <- rep(NA_real_, n)
  since_off <- rep(NA_real_, n)
  on_cnt <- 0; off_cnt <- Inf
  for (i in seq_len(n)) {
    if (state[i]) { on_cnt <- on_cnt + 1; off_cnt <- 0 } else { on_cnt <- 0; off_cnt <- off_cnt + 1 }
    since_on[i] <- on_cnt * step_min / 60
    since_off[i] <- if (is.finite(off_cnt)) off_cnt * step_min / 60 else NA_real_
  }
  since_off[state] <- NA_real_
  data.frame(time = time, since_onset_h = since_on, since_erosion_h = since_off)
}


#' Vertical weighting over the height of the icefall.
#'
#' OPTIONAL - not part of the default correction, see the note in CAP_PARAMS.
#'
#' An icefall is a VERTICAL structure that cuts straight through the sharpest
#' part of the CAP temperature profile. Rauchoecker et al. (2024) measured an
#' inversion depth of ~30 m with ~2 degC/10 m inside it and isothermal
#' stratification above. Assigning one single temperature to a 30-100 m tall
#' icefall therefore overestimates the cold anomaly for the upper part.
#' Use this as a sensitivity experiment ("how much would the result change if
#' the near-surface sharpening were resolved?"), not stacked on top of the
#' valley-scale cap_delta signal by default.
#'
#' Assumed deficit profile above the valley floor:
#'     d(z) = dT_max * max(0, 1 - z / D),  D = CAP_DEPTH_M
#' The factor is the height-average of d(z)/dT_max over [z_base, z_base + H].
#'
#' @param height_m icefall height H
#' @param base_above_floor_m height of the icefall foot above the valley floor
#' @return factor in [0, 1]; 1 would mean the whole fall sits at the floor.
cap_height_factor <- function(height_m = CAP_PARAMS$DEFAULT_ICEFALL_HEIGHT_M,
                              base_above_floor_m = CAP_PARAMS$DEFAULT_ICEFALL_BASE_M,
                              p = CAP_PARAMS) {
  H <- as.numeric(height_m)
  a <- pmax(0, as.numeric(base_above_floor_m))
  D <- p$CAP_DEPTH_M

  H[!is.finite(H) | H <= 0] <- p$DEFAULT_ICEFALL_HEIGHT_M
  a[!is.finite(a)] <- p$DEFAULT_ICEFALL_BASE_M

  b <- pmin(a + H, D)                       # upper integration limit
  inside <- b > a
  # integral of (1 - z/D) from a to b, divided by H
  val <- ifelse(inside, ((b - a) - (b^2 - a^2) / (2 * D)) / H, 0)
  clamp01(val)
}


# ---------------------------------------------------------------------------
# 3) Main function: physical CAP temperature offset [K]
# ---------------------------------------------------------------------------

#' Cold-air-pool temperature offset for an icefall relative to its AWS.
#'
#' The result is NEGATIVE (or zero): the icefall is colder than the
#' lapse-rate-only prediction whenever it sits deeper in the cold pool than
#' its reference station.
#'
#' @param cap_delta_positive pmax(0, icefall_cap_potential - station_cap_potential)
#'   from data/CAP/icefall_station_cap.csv
#' @param inv_grad_C_per_100m,inv_active from the inversion cache
#' @param ff_ms wind speed at the reference station [m/s]
#' @param since_onset_h,since_erosion_h from cap_state_clock()
#' @param height_m,base_above_floor_m icefall geometry (see cap_height_factor)
#' @param cap_pair_confidence from icefall_station_cap.csv
#' @param use_height_weighting FALSE by default (scale mismatch, see CAP_PARAMS);
#'   set TRUE only for the vertical-structure sensitivity run
#' @param K_CAP scaling [K] for cap_delta_positive = 1
#' @return numeric vector of dT_cap in Kelvin (<= 0)
cap_temperature_offset <- function(cap_delta_positive,
                                   inv_grad_C_per_100m = NA_real_,
                                   inv_active = NULL,
                                   ff_ms = NA_real_,
                                   since_onset_h = Inf,
                                   since_erosion_h = NA_real_,
                                   height_m = CAP_PARAMS$DEFAULT_ICEFALL_HEIGHT_M,
                                   base_above_floor_m = CAP_PARAMS$DEFAULT_ICEFALL_BASE_M,
                                   cap_pair_confidence = 1,
                                   use_height_weighting = FALSE,
                                   inv_mode = "hybrid",
                                   K_CAP = CAP_PARAMS$K_CAP,
                                   p = CAP_PARAMS) {

  if (K_CAP < 0 || K_CAP > p$K_CAP_HARD_MAX) {
    stop(sprintf(paste0("K_CAP = %.2f K is outside the physically plausible range ",
                        "[0, %.1f] K. The upper bound is the Gruenloch extreme ",
                        "(Pospichal et al. 2003); ordinary Alpine valleys are ",
                        "%.1f-%.1f K (Rauchoecker et al. 2024)."),
                 K_CAP, p$K_CAP_HARD_MAX, p$K_CAP_MIN_LIT, p$K_CAP_MAX_LIT))
  }
  if (K_CAP > p$K_CAP_MAX_LIT) {
    warning(sprintf(paste0("K_CAP = %.2f K exceeds the range measured in ordinary ",
                           "Alpine valleys (%.1f-%.1f K). Justify this explicitly."),
                    K_CAP, p$K_CAP_MIN_LIT, p$K_CAP_MAX_LIT))
  }

  d <- clamp01(as.numeric(cap_delta_positive)); d[!is.finite(d)] <- 0

  f_inv  <- cap_inversion_factor(inv_grad_C_per_100m, inv_active, mode = inv_mode, p = p)
  f_wind <- cap_wind_factor(ff_ms, p = p)
  f_ramp <- cap_ramp_factor(since_onset_h, since_erosion_h, p = p)
  f_h    <- if (isTRUE(use_height_weighting)) {
    cap_height_factor(height_m, base_above_floor_m, p = p)
  } else 1

  f_conf <- clamp01(as.numeric(cap_pair_confidence)); f_conf[!is.finite(f_conf)] <- 0

  -K_CAP * d * f_inv * f_wind * f_ramp * f_h * f_conf
}


# ---------------------------------------------------------------------------
# 4) Convenience wrapper for the operational model
# ---------------------------------------------------------------------------

#' Add lapse-rate AND CAP-corrected icefall temperature to a model data frame.
#'
#' Replaces the previous one-liner
#'     TL_ice <- TL_station - GAMMA * dz
#'
#' @param df data frame with at least: TL_station, and the columns named in
#'   the arguments below.
#' @param gamma_K_per_m environmental lapse rate (default 0.0065)
#' @return df with new columns dT_lapse, dT_cap, TL_ice_nocap, TL_ice
apply_cap_correction <- function(df,
                                 col_T = "TL_station",
                                 col_dz = "dz_m",
                                 col_cap = "cap_delta_positive",
                                 col_grad = "inv_grad_max_C_per_100m",
                                 col_active = "inv_active",
                                 col_ff = "FF",
                                 col_onset = "since_onset_h",
                                 col_erode = "since_erosion_h",
                                 col_height = "icefall_height_m",
                                 col_base = "icefall_base_above_floor_m",
                                 col_conf = "cap_pair_confidence",
                                 gamma_K_per_m = 0.0065,
                                 K_CAP = CAP_PARAMS$K_CAP,
                                 use_height_weighting = FALSE,
                                 inv_mode = "hybrid") {

  g <- function(nm, default) if (!is.null(df[[nm]])) df[[nm]] else default

  df$dT_lapse <- -gamma_K_per_m * as.numeric(g(col_dz, 0))

  df$dT_cap <- cap_temperature_offset(
    cap_delta_positive  = g(col_cap, 0),
    inv_grad_C_per_100m = g(col_grad, NA_real_),
    inv_active          = g(col_active, NULL),
    ff_ms               = g(col_ff, NA_real_),
    since_onset_h       = g(col_onset, Inf),
    since_erosion_h     = g(col_erode, NA_real_),
    height_m            = g(col_height, CAP_PARAMS$DEFAULT_ICEFALL_HEIGHT_M),
    base_above_floor_m  = g(col_base, CAP_PARAMS$DEFAULT_ICEFALL_BASE_M),
    cap_pair_confidence = g(col_conf, 1),
    use_height_weighting = use_height_weighting,
    inv_mode            = inv_mode,
    K_CAP               = K_CAP
  )

  df$TL_ice_nocap <- as.numeric(df[[col_T]]) + df$dT_lapse
  df$TL_ice       <- df$TL_ice_nocap + df$dT_cap
  df
}


# ---------------------------------------------------------------------------
# 5) Self-test / sanity check  (source this file and call cap_selftest())
# ---------------------------------------------------------------------------

cap_selftest <- function() {
  cat("--- cap_temperature.R self-test ---------------------------------\n")

  cat("\nHeight weighting (CAP depth = ", CAP_PARAMS$CAP_DEPTH_M, " m, base at floor):\n", sep = "")
  for (H in c(5, 10, 20, 30, 50, 100, 150)) {
    cat(sprintf("  H = %3d m -> factor %.3f  (mean offset %.2f K at cap_delta=1)\n",
                H, cap_height_factor(H), CAP_PARAMS$K_CAP * cap_height_factor(H)))
  }

  cat("\nWind gating:\n")
  for (ff in c(0, 0.5, 1.0, 1.5, 2.0, 3.0)) {
    cat(sprintf("  FF = %.1f m/s -> factor %.2f\n", ff, cap_wind_factor(ff)))
  }

  cat("\nBuild-up (tau = ", CAP_PARAMS$TAU_BUILD_H, " h):\n", sep = "")
  for (h in c(0.5, 1, 2, 3, 6, 12)) {
    cat(sprintf("  %4.1f h after onset -> factor %.2f\n", h, cap_ramp_factor(h)))
  }

  cat("\nDEFAULT correction, ideal CAP night (calm 0.3 m/s, strong inversion\n")
  cat("0.6 degC/100 m, 6 h after onset, no height weighting):\n")
  for (d in c(0.142, 0.344, 0.524, 0.781)) {
    dT <- cap_temperature_offset(
      cap_delta_positive = d, inv_grad_C_per_100m = 0.6, inv_active = TRUE,
      ff_ms = 0.3, since_onset_h = 6
    )
    cat(sprintf("  cap_delta_positive = %.3f -> dT_cap = %+.2f K\n", d, dT))
  }
  cat("\nSAME with the optional height weighting (30 m fall on the floor):\n")
  for (d in c(0.142, 0.344, 0.524, 0.781)) {
    dT <- cap_temperature_offset(
      cap_delta_positive = d, inv_grad_C_per_100m = 0.6, inv_active = TRUE,
      ff_ms = 0.3, since_onset_h = 6, height_m = 30,
      use_height_weighting = TRUE
    )
    cat(sprintf("  cap_delta_positive = %.3f -> dT_cap = %+.2f K\n", d, dT))
  }
  cat("\n(0.142 / 0.344 / 0.524 / 0.781 = median / p75 / p90 / max of the 303\n",
      " icefall-station pairs in data/CAP/icefall_station_cap.csv)\n", sep = "")
  invisible(TRUE)
}
