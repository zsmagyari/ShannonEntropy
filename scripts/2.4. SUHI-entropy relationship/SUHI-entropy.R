library(terra)
library(mgcv)

locs <- c("BRASOV","CLUJ_NAPOCA","CONSTANTA","CRAIOVA","IASI","TIMISOARA")

interval <- "2021-2025"

BASE_SCREEN_DIR <- file.path("ANALYSIS", "sensitivity")
BASE_DW_DIR     <- "DW"
BASE_LST_DIR    <- "LST"
BASE_DEM_DIR    <- "DEMs"
BASE_ALIGN_DIR  <- file.path("ANALYSIS", "aligned")

BASE_OUT_DIR    <- file.path("ANALYSIS", "final_screened")

TMP_DIR <- "F:/RTMP"
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)
terraOptions(tempdir = TMP_DIR, todisk = TRUE, memfrac = 0.60, progress = 1)

TH_CORE   <- 0.6
TH_WATER  <- 0.5
TH_SETTLE <- 0.6
OTHER_SETT_BUF_M <- 500
DIRECTIONS <- 8
KEEP_ONLY_CONNECTED_TO_CORE <- TRUE
USE_URBAN_PIXELS_THRESHOLD  <- TRUE

REF_RURAL <- list(
  TH_RURAL    = 0.1,
  K_RURAL_MAX = 2.0,
  ELEV_TOL_M  = 1000
)

GRID_M <- 300
JITTER_WITHIN_CELL <- TRUE
JITTER_SEED_FIXED  <- 1
MIN_URBAN_PX <- 200
MIN_RURAL_PX <- 200
MIN_POINTS   <- 80

K_ENTROPY <- 8
KXY_MANUAL <- c(
  BRASOV      = 140,
  BUCURESTI   = 300,
  CLUJ_NAPOCA = 180,
  CONSTANTA   = 150,
  CRAIOVA     = 150,
  IASI        = 250,
  TIMISOARA   = 300
)

choose_kxy <- function(n_points, k_min = 80, k_max = 250) {
  k <- round(5 * sqrt(n_points))
  k <- max(k_min, k)
  k <- min(k_max, k)
  k
}

ZERO_TOL <- 1e-8
WINDOW_SELECTION_MODE <- "all_tested"
WINDOW_COMPETITIVE_DELTA <- 2

find_one <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(f) == 0) stop("Missing file. Pattern: ", pattern, " | dir: ", dir)
  if (length(f) > 1) message("Multiple matches for pattern ", pattern, " -> using first: ", basename(f[1]))
  f[1]
}

make_pts_grid <- function(template_r, grid_m = 300, jitter = TRUE, seed = 1) {
  gridr <- rast(ext(template_r), resolution = grid_m, crs = crs(template_r))
  pts0  <- as.points(gridr, values = FALSE)
  if (!jitter) return(pts0)

  xy <- crds(pts0, df = TRUE)
  set.seed(seed)
  xy$x <- xy$x + runif(nrow(xy), -grid_m / 2, grid_m / 2)
  xy$y <- xy$y + runif(nrow(xy), -grid_m / 2, grid_m / 2)

  e <- ext(template_r)
  xy$x <- pmin(pmax(xy$x, e[1]), e[2])
  xy$y <- pmin(pmax(xy$y, e[3]), e[4])
  vect(xy, geom = c("x", "y"), crs = crs(template_r))
}

g_median <- function(r) {
  out <- try(global(r, fun = "median", na.rm = TRUE), silent = TRUE)
  if (!inherits(out, "try-error") && is.finite(out[1, 1])) return(as.numeric(out[1, 1]))
  q <- global(r, fun = function(x, ...) stats::quantile(x, probs = 0.5, na.rm = TRUE), na.rm = TRUE)
  as.numeric(q[1, 1])
}

count_cells <- function(mask01) {
  v <- global(mask01, "sum", na.rm = TRUE)[1, 1]
  if (!is.finite(v)) 0 else as.numeric(v)
}

aicc_from_model <- function(model) {
  aic <- AIC(model)
  ll  <- logLik(model)
  k   <- attr(ll, "df")
  n   <- stats::nobs(model)
  if (!is.finite(aic) || !is.finite(k) || !is.finite(n)) return(NA_real_)
  if ((n - k - 1) <= 0) return(NA_real_)
  as.numeric(aic + (2 * k * (k + 1)) / (n - k - 1))
}

safe_s_table <- function(gam_summary, row_idx, col_name) {
  tryCatch(as.numeric(gam_summary$s.table[row_idx, col_name]), error = function(e) NA_real_)
}

safe_edf <- function(gam_summary, row_idx) {
  nm <- colnames(gam_summary$s.table)
  if ("edf" %in% nm) return(safe_s_table(gam_summary, row_idx, "edf"))
  NA_real_
}

build_masks <- function(b_lst, w_lst, dem_lst,
                        TH_CORE, TH_WATER,
                        TH_URBAN, BUFFER_M,
                        TH_RURAL, K_RURAL_MAX, ELEV_TOL_M,
                        TH_SETTLE, OTHER_SETT_BUF_M,
                        DIRECTIONS,
                        KEEP_ONLY_CONNECTED_TO_CORE = TRUE,
                        USE_URBAN_PIXELS_THRESHOLD = TRUE) {

  core_cand <- ifel((b_lst >= TH_CORE) & (w_lst < TH_WATER), 1, NA)
  p_core <- patches(core_cand, directions = DIRECTIONS)
  fr <- as.data.frame(freq(p_core))
  fr <- fr[!is.na(fr$value) & fr$value != 0, , drop = FALSE]
  if (nrow(fr) == 0) stop("No core candidate pixels found.")

  core_id <- fr$value[which.max(fr$count)]
  core <- ifel(p_core == core_id, 1, NA)
  names(core) <- "urban_core"

  d <- distance(core)

  core_cells <- global(core, "sum", na.rm = TRUE)[1, 1]
  cell_area_m2 <- abs(res(b_lst)[1] * res(b_lst)[2])
  core_area_m2 <- core_cells * cell_area_m2
  r_eq_m <- sqrt(core_area_m2 / pi)

  core_dem <- mask(dem_lst, core)
  core_dem_mean <- global(core_dem, "mean", na.rm = TRUE)[1, 1]
  if (!is.finite(core_dem_mean)) stop("Core mean elevation is NA.")

  if (USE_URBAN_PIXELS_THRESHOLD) {
    urban_pix <- ifel((d <= BUFFER_M) & (b_lst >= TH_URBAN) & (w_lst < TH_WATER), 1, NA)
  } else {
    urban_pix <- ifel((d <= BUFFER_M) & (w_lst < TH_WATER), 1, NA)
  }
  names(urban_pix) <- "urban_pix_raw"

  if (KEEP_ONLY_CONNECTED_TO_CORE && USE_URBAN_PIXELS_THRESHOLD) {
    up <- patches(ifel(!is.na(urban_pix), 1, NA), directions = DIRECTIONS)
    ids <- unique(na.omit(values(mask(up, core))))
    if (length(ids) == 0) {
      urban_final <- urban_pix
    } else {
      urban_final <- ifel(up %in% ids, 1, NA)
    }
  } else {
    urban_final <- urban_pix
  }

  urban_final <- cover(urban_final, core)
  names(urban_final) <- "urban"

  settle_cand <- ifel((b_lst >= TH_SETTLE) & (w_lst < TH_WATER), 1, NA)
  sp <- patches(settle_cand, directions = DIRECTIONS)
  sfr <- as.data.frame(freq(sp))
  sfr <- sfr[!is.na(sfr$value) & sfr$value != 0, , drop = FALSE]

  target_ids <- unique(na.omit(values(mask(sp, core))))
  other_ids  <- setdiff(sfr$value, target_ids)

  other_sett <- ifel(sp %in% other_ids, 1, NA)
  has_other <- is.finite(global(other_sett, "max", na.rm = TRUE)[1, 1])

  if (has_other) {
    dist_other <- distance(other_sett)
    keep_outside_other <- ifel(dist_other > OTHER_SETT_BUF_M, 1, NA)
  } else {
    keep_outside_other <- ifel(!is.na(b_lst), 1, NA)
  }

  RURAL_MAX_M <- as.numeric(K_RURAL_MAX) * r_eq_m

  rural <- ifel(
    (d > BUFFER_M) &
      (d <= RURAL_MAX_M) &
      (b_lst <= TH_RURAL) &
      (w_lst < TH_WATER) &
      (abs(dem_lst - core_dem_mean) <= ELEV_TOL_M),
    1, NA
  )
  rural <- mask(rural, keep_outside_other)
  names(rural) <- "rural"

  list(
    core = core,
    urban = urban_final,
    rural = rural
  )
}

read_screening_long <- function(city_key, interval, base_dir) {
  pat <- sprintf("^RO_%s_.*results_LONG_%s\\.csv$", city_key, interval)
  f <- list.files(base_dir, pattern = pat, full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  if (length(f) == 0) stop("Missing screening LONG CSV for city: ", city_key)
  if (length(f) > 1) message("Multiple screening files for ", city_key, " -> using first: ", basename(f[1]))
  read.csv(f[1], stringsAsFactors = FALSE)
}

screen_candidates <- function(screen_df, zero_tol = 1e-8,
                              window_selection_mode = "all_tested",
                              window_competitive_delta = 2) {
  required_cols <- c("city", "window", "TH_URBAN", "TH_RURAL", "BUFFER_M", "K_RURAL_MAX", "ELEV_TOL_M", "best_aicc")
  miss <- setdiff(required_cols, names(screen_df))
  if (length(miss) > 0) stop("Missing screening columns: ", paste(miss, collapse = ", "))

  min_aicc <- min(screen_df$best_aicc, na.rm = TRUE)
  screen_df$delta_to_city_min <- screen_df$best_aicc - min_aicc
  zero_set <- screen_df[is.finite(screen_df$delta_to_city_min) & abs(screen_df$delta_to_city_min) <= zero_tol, , drop = FALSE]
  if (nrow(zero_set) == 0) stop("No rows in zero-ΔAICc set.")

  urban_pairs <- unique(zero_set[, c("TH_URBAN", "BUFFER_M")])
  if (nrow(urban_pairs) != 1) stop("Zero-ΔAICc set does not yield a unique urban-side setting.")

  rural_set <- unique(zero_set[, c("TH_RURAL", "K_RURAL_MAX", "ELEV_TOL_M")])
  rural_set <- rural_set[order(rural_set$TH_RURAL, rural_set$K_RURAL_MAX, rural_set$ELEV_TOL_M), , drop = FALSE]

  window_summary <- aggregate(best_aicc ~ window, data = screen_df, FUN = function(x) min(x, na.rm = TRUE))
  names(window_summary)[names(window_summary) == "best_aicc"] <- "window_min_best_aicc"
  window_summary$window_m <- as.numeric(gsub("[^0-9]", "", window_summary$window))
  window_summary$delta_window_min_to_city_min <- window_summary$window_min_best_aicc - min_aicc

  zero_n_by_window <- aggregate(rep(1, nrow(zero_set)) ~ window, data = zero_set, FUN = sum)
  names(zero_n_by_window)[2] <- "n_zero_rows"
  window_summary <- merge(window_summary, zero_n_by_window, by = "window", all.x = TRUE, sort = FALSE)
  window_summary$n_zero_rows[is.na(window_summary$n_zero_rows)] <- 0L

  fixed_urban_rows <- screen_df[
    screen_df$TH_URBAN == urban_pairs$TH_URBAN[1] &
      screen_df$BUFFER_M == urban_pairs$BUFFER_M[1],
    , drop = FALSE
  ]
  fixed_urban_min_by_window <- aggregate(best_aicc ~ window, data = fixed_urban_rows, FUN = function(x) min(x, na.rm = TRUE))
  names(fixed_urban_min_by_window)[2] <- "fixed_urban_window_min_best_aicc"
  window_summary <- merge(window_summary, fixed_urban_min_by_window, by = "window", all.x = TRUE, sort = FALSE)
  window_summary$delta_fixed_urban_window_min_to_city_min <- window_summary$fixed_urban_window_min_best_aicc - min_aicc
  window_summary <- window_summary[order(window_summary$window_m), , drop = FALSE]

  all_windows <- sort(unique(window_summary$window_m))
  if (length(all_windows) == 0) stop("No windows parsed from screening file.")

  if (window_selection_mode == "all_tested") {
    selected_windows <- all_windows
  } else if (window_selection_mode == "competitive") {
    selected_windows <- sort(unique(window_summary$window_m[window_summary$delta_window_min_to_city_min <= window_competitive_delta]))
    if (length(selected_windows) == 0) selected_windows <- all_windows
  } else {
    stop("Unsupported WINDOW_SELECTION_MODE.")
  }

  list(
    min_aicc = min_aicc,
    zero_set = zero_set,
    urban = urban_pairs,
    rural_set = rural_set,
    window_summary = window_summary,
    selected_windows = selected_windows
  )
}

pick_reference_rural <- function(rural_set, ref_rural) {
  hit <- which(
    rural_set$TH_RURAL == ref_rural$TH_RURAL &
      rural_set$K_RURAL_MAX == ref_rural$K_RURAL_MAX &
      rural_set$ELEV_TOL_M == ref_rural$ELEV_TOL_M
  )
  if (length(hit) > 0) return(rural_set[hit[1], , drop = FALSE])

  ord <- order(rural_set$TH_RURAL,
               abs(rural_set$K_RURAL_MAX - median(rural_set$K_RURAL_MAX)),
               rural_set$ELEV_TOL_M)
  rural_set[ord[1], , drop = FALSE]
}

prepare_city_inputs <- function(CITY_KEY, interval,
                                BASE_DW_DIR, BASE_LST_DIR, BASE_DEM_DIR, BASE_ALIGN_DIR) {

  DW_DIR   <- file.path(BASE_DW_DIR, CITY_KEY)
  LST_DIR  <- file.path(BASE_LST_DIR, CITY_KEY)
  DEM_PATH <- file.path(BASE_DEM_DIR, sprintf("%s_dem.tif", CITY_KEY))
  ALIGNED_DIR <- file.path(BASE_ALIGN_DIR, CITY_KEY)

  built_path <- find_one(DW_DIR, sprintf("^RO_%s_DW_JJA_MULTIYEAR_%s_.*_built\\.tif$", CITY_KEY, interval))
  water_path <- find_one(DW_DIR, sprintf("^RO_%s_DW_JJA_MULTIYEAR_%s_.*_water\\.tif$", CITY_KEY, interval))
  lst_path   <- find_one(LST_DIR, sprintf("(?i)^RO_%s_.*LST_HARMONIZED.*JJA.*MULTIYEAR.*%s.*\\.tif$", CITY_KEY, interval))

  b <- rast(built_path)
  w <- rast(water_path)
  if (!compareGeom(b, w, stopOnError = FALSE)) w <- resample(w, b, method = "bilinear")

  if (!file.exists(DEM_PATH)) stop("Missing DEM file: ", DEM_PATH)
  dem <- rast(DEM_PATH)[[1]]
  if (!same.crs(dem, b)) dem <- project(dem, b, method = "bilinear")
  if (!compareGeom(dem, b, stopOnError = FALSE)) dem <- resample(dem, b, method = "bilinear")

  lst <- rast(lst_path)
  if (!same.crs(b, lst) || !compareGeom(b, lst, stopOnError = FALSE)) {
    b_lst   <- project(b,   lst, method = "bilinear")
    w_lst   <- project(w,   lst, method = "bilinear")
    dem_lst <- project(dem, lst, method = "bilinear")
  } else {
    b_lst <- b
    w_lst <- w
    dem_lst <- dem
  }

  list(
    b_lst = b_lst,
    w_lst = w_lst,
    dem_lst = dem_lst,
    lst = lst,
    aligned_dir = ALIGNED_DIR
  )
}

load_entropy_raster <- function(CITY_KEY, interval, window_m, aligned_dir, lst_template) {
  ent_pat1 <- sprintf("(?i)^RO_%s_ENTROPY_%s_W%dm_ALIGNED_TO_SUHI\\.tif$", CITY_KEY, interval, window_m)
  ent_path <- list.files(aligned_dir, pattern = ent_pat1, full.names = TRUE)
  if (length(ent_path) == 0) {
    ent_pat2 <- sprintf("(?i)%s.*W%dm.*ALIGNED.*\\.tif$", interval, window_m)
    ent_path <- list.files(aligned_dir, pattern = ent_pat2, full.names = TRUE)
  }
  if (length(ent_path) == 0) stop("Missing entropy aligned file for W", window_m, "m in ", aligned_dir)
  ent <- rast(ent_path[1])
  if (!same.crs(ent, lst_template) || !compareGeom(ent, lst_template, stopOnError = FALSE)) {
    ent <- project(ent, lst_template, method = "bilinear")
  }
  names(ent) <- "entropy"
  ent
}

run_one_analysis <- function(CITY_KEY, interval, window_m,
                             city_inputs,
                             TH_URBAN, BUFFER_M,
                             TH_RURAL, K_RURAL_MAX, ELEV_TOL_M,
                             TH_CORE, TH_WATER, TH_SETTLE, OTHER_SETT_BUF_M,
                             DIRECTIONS,
                             KEEP_ONLY_CONNECTED_TO_CORE,
                             USE_URBAN_PIXELS_THRESHOLD,
                             GRID_M, JITTER_WITHIN_CELL, JITTER_SEED_FIXED,
                             MIN_URBAN_PX, MIN_RURAL_PX, MIN_POINTS,
                             K_ENTROPY, KXY_MANUAL,
                             pts_override = NULL,
                             common_ok_override = NULL) {

  b_lst <- city_inputs$b_lst
  w_lst <- city_inputs$w_lst
  dem_lst <- city_inputs$dem_lst
  lst <- city_inputs$lst
  ent <- load_entropy_raster(CITY_KEY, interval, window_m, city_inputs$aligned_dir, lst)

  masks <- build_masks(
    b_lst = b_lst, w_lst = w_lst, dem_lst = dem_lst,
    TH_CORE = TH_CORE, TH_WATER = TH_WATER,
    TH_URBAN = TH_URBAN, BUFFER_M = BUFFER_M,
    TH_RURAL = TH_RURAL, K_RURAL_MAX = K_RURAL_MAX, ELEV_TOL_M = ELEV_TOL_M,
    TH_SETTLE = TH_SETTLE, OTHER_SETT_BUF_M = OTHER_SETT_BUF_M,
    DIRECTIONS = DIRECTIONS,
    KEEP_ONLY_CONNECTED_TO_CORE = KEEP_ONLY_CONNECTED_TO_CORE,
    USE_URBAN_PIXELS_THRESHOLD = USE_URBAN_PIXELS_THRESHOLD
  )

  urban <- masks$urban
  rural <- masks$rural

  n_urban_px <- count_cells(urban)
  n_rural_px <- count_cells(rural)
  if (n_urban_px < MIN_URBAN_PX) stop("Urban mask too small: ", n_urban_px)
  if (n_rural_px < MIN_RURAL_PX) stop("Rural mask too small: ", n_rural_px)

  med_rural <- g_median(mask(lst, rural))
  suhi <- lst - med_rural
  med_urban_lst <- g_median(mask(lst, urban))
  suhi_delta_median <- med_urban_lst - med_rural

  pts_all <- pts_override
  if (is.null(pts_all)) {
    pts_all <- make_pts_grid(lst, grid_m = GRID_M, jitter = JITTER_WITHIN_CELL, seed = JITTER_SEED_FIXED)
  }

  u_at_pts <- extract(urban, pts_all)[, 2]
  idx_u <- which(is.finite(u_at_pts) & u_at_pts > 0)
  if (length(idx_u) < MIN_POINTS) stop("Too few urban points after thinning: ", length(idx_u))
  pts_u <- pts_all[idx_u]
  xy_u <- crds(pts_u, df = TRUE)

  suhi_u <- extract(suhi, pts_u)[, 2]
  ent_u  <- extract(ent,  pts_u)[, 2]

  ok <- is.finite(suhi_u) & is.finite(ent_u)
  if (!is.null(common_ok_override)) {
    if (length(common_ok_override) != length(ok)) stop("common_ok_override length does not match extracted points.")
    ok <- ok & common_ok_override
  }

  n_pts <- sum(ok)
  if (n_pts < MIN_POINTS) stop("Too few valid points: ", n_pts)

  df <- data.frame(
    x = xy_u$x[ok],
    y = xy_u$y[ok],
    suhi = suhi_u[ok],
    entropy = ent_u[ok]
  )

  rho <- suppressWarnings(cor(df$suhi, df$entropy, method = "spearman"))
  r   <- suppressWarnings(cor(df$suhi, df$entropy, method = "pearson"))

  entropy_iqr   <- IQR(df$entropy, na.rm = TRUE)
  entropy_sd    <- sd(df$entropy, na.rm = TRUE)
  entropy_range <- diff(range(df$entropy, na.rm = TRUE))
  suhi_iqr      <- IQR(df$suhi, na.rm = TRUE)
  suhi_sd       <- sd(df$suhi, na.rm = TRUE)
  suhi_range    <- diff(range(df$suhi, na.rm = TRUE))

  gam1 <- mgcv::gam(suhi ~ s(entropy, k = K_ENTROPY), data = df, method = "REML")

  k_xy_auto <- choose_kxy(nrow(df))
  k_xy <- k_xy_auto
  k_xy_mode <- "adaptive"
  if (CITY_KEY %in% names(KXY_MANUAL) && is.finite(KXY_MANUAL[[CITY_KEY]])) {
    k_xy <- as.integer(KXY_MANUAL[[CITY_KEY]])
    k_xy_mode <- "manual"
  }

  gam2 <- mgcv::gam(
    suhi ~ s(entropy, k = K_ENTROPY) + s(x, y, k = k_xy),
    data = df,
    method = "REML"
  )

  s1 <- summary(gam1)
  s2 <- summary(gam2)
  gam1_aic  <- AIC(gam1)
  gam2_aic  <- AIC(gam2)
  gam1_aicc <- aicc_from_model(gam1)
  gam2_aicc <- aicc_from_model(gam2)

  data.frame(
    city = CITY_KEY,
    interval = interval,
    window_m = window_m,
    TH_CORE = TH_CORE,
    TH_WATER = TH_WATER,
    TH_URBAN = TH_URBAN,
    BUFFER_M = BUFFER_M,
    TH_RURAL = TH_RURAL,
    K_RURAL_MAX = K_RURAL_MAX,
    ELEV_TOL_M = ELEV_TOL_M,
    TH_SETTLE = TH_SETTLE,
    OTHER_SETT_BUF_M = OTHER_SETT_BUF_M,
    n_urban_px = n_urban_px,
    n_rural_px = n_rural_px,
    n_points = n_pts,
    lst_urban_median = med_urban_lst,
    lst_rural_median = med_rural,
    suhi_delta_median_C = suhi_delta_median,
    spearman_rho = rho,
    pearson_r = r,
    entropy_iqr = entropy_iqr,
    entropy_sd = entropy_sd,
    entropy_range = entropy_range,
    suhi_iqr = suhi_iqr,
    suhi_sd = suhi_sd,
    suhi_range = suhi_range,
    k_xy_used = k_xy,
    k_xy_mode = k_xy_mode,
    k_xy_auto = k_xy_auto,
    gam1_aic = gam1_aic,
    gam1_aicc = gam1_aicc,
    gam1_dev_expl = s1$dev.expl,
    gam1_r2_adj = s1$r.sq,
    gam1_edf_entropy = safe_edf(s1, 1),
    gam1_p_smooth_entropy = safe_s_table(s1, 1, "p-value"),
    gam2_aic = gam2_aic,
    gam2_aicc = gam2_aicc,
    gam2_dev_expl = s2$dev.expl,
    gam2_r2_adj = s2$r.sq,
    gam2_edf_entropy = safe_edf(s2, 1),
    gam2_edf_xy = safe_edf(s2, 2),
    gam2_p_smooth_entropy = safe_s_table(s2, 1, "p-value"),
    gam2_p_smooth_xy = safe_s_table(s2, 2, "p-value"),
    row.names = NULL
  )
}

dir.create(BASE_OUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (CITY_KEY in locs) {
  OUT_DIR <- file.path(BASE_OUT_DIR, CITY_KEY)
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

  screen_df <- read_screening_long(CITY_KEY, interval, BASE_SCREEN_DIR)
  cand <- screen_candidates(
    screen_df,
    zero_tol = ZERO_TOL,
    window_selection_mode = WINDOW_SELECTION_MODE,
    window_competitive_delta = WINDOW_COMPETITIVE_DELTA
  )

  zero_set <- cand$zero_set
  urban_set <- cand$urban
  rural_set <- cand$rural_set
  window_summary <- cand$window_summary
  ref_rural <- pick_reference_rural(rural_set, REF_RURAL)

  TH_URBAN_MAIN <- urban_set$TH_URBAN[1]
  BUFFER_MAIN   <- urban_set$BUFFER_M[1]
  candidate_windows <- cand$selected_windows
  if (length(candidate_windows) == 0) stop("No candidate windows parsed from screening stage for ", CITY_KEY)

  zero_set$window_m <- as.numeric(gsub("[^0-9]", "", zero_set$window))
  write.csv(
    zero_set,
    file.path(OUT_DIR, sprintf("RO_%s_screening_zero_set_%s.csv", CITY_KEY, interval)),
    row.names = FALSE
  )

  window_summary$city <- CITY_KEY
  window_summary$interval <- interval
  write.csv(
    window_summary,
    file.path(OUT_DIR, sprintf("RO_%s_screening_window_summary_%s.csv", CITY_KEY, interval)),
    row.names = FALSE
  )

  city_inputs <- prepare_city_inputs(CITY_KEY, interval, BASE_DW_DIR, BASE_LST_DIR, BASE_DEM_DIR, BASE_ALIGN_DIR)

  pts_all_city <- make_pts_grid(city_inputs$lst, grid_m = GRID_M, jitter = JITTER_WITHIN_CELL, seed = JITTER_SEED_FIXED)

  masks_ref <- build_masks(
    b_lst = city_inputs$b_lst, w_lst = city_inputs$w_lst, dem_lst = city_inputs$dem_lst,
    TH_CORE = TH_CORE, TH_WATER = TH_WATER,
    TH_URBAN = TH_URBAN_MAIN, BUFFER_M = BUFFER_MAIN,
    TH_RURAL = ref_rural$TH_RURAL[1], K_RURAL_MAX = ref_rural$K_RURAL_MAX[1], ELEV_TOL_M = ref_rural$ELEV_TOL_M[1],
    TH_SETTLE = TH_SETTLE, OTHER_SETT_BUF_M = OTHER_SETT_BUF_M,
    DIRECTIONS = DIRECTIONS,
    KEEP_ONLY_CONNECTED_TO_CORE = KEEP_ONLY_CONNECTED_TO_CORE,
    USE_URBAN_PIXELS_THRESHOLD = USE_URBAN_PIXELS_THRESHOLD
  )

  urban_ref <- masks_ref$urban
  lst_ref <- city_inputs$lst
  med_rural_ref <- g_median(mask(lst_ref, masks_ref$rural))
  suhi_ref <- lst_ref - med_rural_ref

  u_at_pts <- extract(urban_ref, pts_all_city)[, 2]
  idx_u <- which(is.finite(u_at_pts) & u_at_pts > 0)
  if (length(idx_u) < MIN_POINTS) stop("Too few urban points after thinning for main comparison in ", CITY_KEY)
  pts_u_city <- pts_all_city[idx_u]

  suhi_u_city <- extract(suhi_ref, pts_u_city)[, 2]
  common_ok <- is.finite(suhi_u_city)
  for (w_m in candidate_windows) {
    ent_tmp <- load_entropy_raster(CITY_KEY, interval, w_m, city_inputs$aligned_dir, city_inputs$lst)
    ent_u_tmp <- extract(ent_tmp, pts_u_city)[, 2]
    common_ok <- common_ok & is.finite(ent_u_tmp)
  }
  if (sum(common_ok) < MIN_POINTS) stop("Common point set across candidate windows too small in ", CITY_KEY)

  main_results <- list()
  for (w_m in candidate_windows) {
    res_main <- run_one_analysis(
      CITY_KEY = CITY_KEY, interval = interval, window_m = w_m,
      city_inputs = city_inputs,
      TH_URBAN = TH_URBAN_MAIN, BUFFER_M = BUFFER_MAIN,
      TH_RURAL = ref_rural$TH_RURAL[1], K_RURAL_MAX = ref_rural$K_RURAL_MAX[1], ELEV_TOL_M = ref_rural$ELEV_TOL_M[1],
      TH_CORE = TH_CORE, TH_WATER = TH_WATER, TH_SETTLE = TH_SETTLE, OTHER_SETT_BUF_M = OTHER_SETT_BUF_M,
      DIRECTIONS = DIRECTIONS,
      KEEP_ONLY_CONNECTED_TO_CORE = KEEP_ONLY_CONNECTED_TO_CORE,
      USE_URBAN_PIXELS_THRESHOLD = USE_URBAN_PIXELS_THRESHOLD,
      GRID_M = GRID_M, JITTER_WITHIN_CELL = JITTER_WITHIN_CELL, JITTER_SEED_FIXED = JITTER_SEED_FIXED,
      MIN_URBAN_PX = MIN_URBAN_PX, MIN_RURAL_PX = MIN_RURAL_PX, MIN_POINTS = MIN_POINTS,
      K_ENTROPY = K_ENTROPY, KXY_MANUAL = KXY_MANUAL,
      pts_override = pts_all_city,
      common_ok_override = common_ok
    )
    res_main$analysis_role <- "MAIN_REFERENCE"
    main_results[[length(main_results) + 1]] <- res_main
  }

  main_summary_df <- do.call(rbind, main_results)
  main_summary_df$delta_gam2_aicc_city <- main_summary_df$gam2_aicc - min(main_summary_df$gam2_aicc, na.rm = TRUE)
  main_summary_df$delta_gam1_aicc_city <- main_summary_df$gam1_aicc - min(main_summary_df$gam1_aicc, na.rm = TRUE)
  main_summary_df$rank_gam2_aicc_city  <- rank(main_summary_df$gam2_aicc, ties.method = "min")
  main_summary_df$rank_abs_rho_city    <- rank(-abs(main_summary_df$spearman_rho), ties.method = "min")
  main_summary_df$rank_gam2_devexpl_city <- rank(-main_summary_df$gam2_dev_expl, ties.method = "min")
  main_summary_df$is_gam2_best_city <- main_summary_df$rank_gam2_aicc_city == 1

  write.csv(
    main_summary_df,
    file.path(OUT_DIR, sprintf("RO_%s_main_GAM_summary_%s.csv", CITY_KEY, interval)),
    row.names = FALSE
  )

  rob_results <- list()
  for (w_m in candidate_windows) {
    u_at_pts_w <- extract(urban_ref, pts_all_city)[, 2]
    idx_u_w <- which(is.finite(u_at_pts_w) & u_at_pts_w > 0)
    pts_u_w <- pts_all_city[idx_u_w]
    ent_w <- load_entropy_raster(CITY_KEY, interval, w_m, city_inputs$aligned_dir, city_inputs$lst)
    ent_u_w <- extract(ent_w, pts_u_w)[, 2]
    common_ok_w <- is.finite(ent_u_w)

    for (ii in seq_len(nrow(rural_set))) {
      rr <- rural_set[ii, , drop = FALSE]
      masks_rr <- build_masks(
        b_lst = city_inputs$b_lst, w_lst = city_inputs$w_lst, dem_lst = city_inputs$dem_lst,
        TH_CORE = TH_CORE, TH_WATER = TH_WATER,
        TH_URBAN = TH_URBAN_MAIN, BUFFER_M = BUFFER_MAIN,
        TH_RURAL = rr$TH_RURAL[1], K_RURAL_MAX = rr$K_RURAL_MAX[1], ELEV_TOL_M = rr$ELEV_TOL_M[1],
        TH_SETTLE = TH_SETTLE, OTHER_SETT_BUF_M = OTHER_SETT_BUF_M,
        DIRECTIONS = DIRECTIONS,
        KEEP_ONLY_CONNECTED_TO_CORE = KEEP_ONLY_CONNECTED_TO_CORE,
        USE_URBAN_PIXELS_THRESHOLD = USE_URBAN_PIXELS_THRESHOLD
      )
      med_rural_rr <- g_median(mask(city_inputs$lst, masks_rr$rural))
      suhi_rr <- city_inputs$lst - med_rural_rr
      suhi_u_rr <- extract(suhi_rr, pts_u_w)[, 2]
      common_ok_w <- common_ok_w & is.finite(suhi_u_rr)
    }

    if (sum(common_ok_w) < MIN_POINTS) stop("Common robustness point set too small in ", CITY_KEY, " | W=", w_m)

    for (ii in seq_len(nrow(rural_set))) {
      rr <- rural_set[ii, , drop = FALSE]
      is_ref <- rr$TH_RURAL == ref_rural$TH_RURAL[1] &&
        rr$K_RURAL_MAX == ref_rural$K_RURAL_MAX[1] &&
        rr$ELEV_TOL_M == ref_rural$ELEV_TOL_M[1]
      if (is_ref) next

      res_rob <- run_one_analysis(
        CITY_KEY = CITY_KEY, interval = interval, window_m = w_m,
        city_inputs = city_inputs,
        TH_URBAN = TH_URBAN_MAIN, BUFFER_M = BUFFER_MAIN,
        TH_RURAL = rr$TH_RURAL[1], K_RURAL_MAX = rr$K_RURAL_MAX[1], ELEV_TOL_M = rr$ELEV_TOL_M[1],
        TH_CORE = TH_CORE, TH_WATER = TH_WATER, TH_SETTLE = TH_SETTLE, OTHER_SETT_BUF_M = OTHER_SETT_BUF_M,
        DIRECTIONS = DIRECTIONS,
        KEEP_ONLY_CONNECTED_TO_CORE = KEEP_ONLY_CONNECTED_TO_CORE,
        USE_URBAN_PIXELS_THRESHOLD = USE_URBAN_PIXELS_THRESHOLD,
        GRID_M = GRID_M, JITTER_WITHIN_CELL = JITTER_WITHIN_CELL, JITTER_SEED_FIXED = JITTER_SEED_FIXED,
        MIN_URBAN_PX = MIN_URBAN_PX, MIN_RURAL_PX = MIN_RURAL_PX, MIN_POINTS = MIN_POINTS,
        K_ENTROPY = K_ENTROPY, KXY_MANUAL = KXY_MANUAL,
        pts_override = pts_all_city,
        common_ok_override = common_ok_w
      )
      res_rob$analysis_role <- "ROBUSTNESS_RURAL"
      rob_results[[length(rob_results) + 1]] <- res_rob
    }
  }

  if (length(rob_results) > 0) {
    rob_summary_df <- do.call(rbind, rob_results)
    write.csv(
      rob_summary_df,
      file.path(OUT_DIR, sprintf("RO_%s_robustness_GAM_summary_%s.csv", CITY_KEY, interval)),
      row.names = FALSE
    )
  }
}