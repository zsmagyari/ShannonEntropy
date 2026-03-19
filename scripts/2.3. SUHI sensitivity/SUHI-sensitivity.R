library(terra)

locs <- c("BRASOV","BUCURESTI","CLUJ_NAPOCA","CONSTANTA","CRAIOVA","IASI","TIMISOARA")

for (CITY_KEY in locs)
{
  interval <- "2021-2025"
  
  DW_DIR   <- file.path("DW", CITY_KEY)
  DEM_PATH <- file.path("DEMs", sprintf("%s_dem.tif", CITY_KEY))
  LST_DIR  <- file.path("LST", CITY_KEY)
  
  ALIGNED_DIR <- file.path("ANALYSIS", "aligned", CITY_KEY)
  
  OUT_DIR <- file.path("ANALYSIS", "sensitivity", CITY_KEY)
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  
  tmp_dir <- "F:/RTMP"
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  terraOptions(tempdir = tmp_dir, todisk = TRUE, memfrac = 0.60, progress = 1)
  
  TH_CORE   <- 0.6
  TH_WATER  <- 0.5
  TH_SETTLE <- 0.6
  OTHER_SETT_BUF_M <- 500
  DIRECTIONS <- 8
  KEEP_ONLY_CONNECTED_TO_CORE <- TRUE
  USE_URBAN_PIXELS_THRESHOLD  <- TRUE
  
  WINDOWS_M <- c(90, 150, 300, 600)
  
  GRID_M <- 300
  JITTER_WITHIN_CELL <- TRUE
  JITTER_SEED_FIXED  <- 1
  
  SENS_GRID <- expand.grid(
    TH_URBAN    = c(0.2, 0.3, 0.4),
    TH_RURAL    = c(0.1, 0.2),
    BUFFER_M    = c(500, 1000, 1500),
    K_RURAL_MAX = c(1.5, 2, 3),
    ELEV_TOL_M  = c(1000, 2000),
    stringsAsFactors = FALSE
  )
  
  find_one <- function(dir, pattern) {
    f <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
    if (length(f) == 0) stop("Missing file. Pattern: ", pattern, " | dir: ", dir)
    if (length(f) > 1) message("Multiple matches for pattern ", pattern, " -> using first: ", basename(f[1]))
    f[1]
  }
  
  built_path <- find_one(DW_DIR, sprintf("^RO_%s_DW_JJA_MULTIYEAR_%s_.*_built\\.tif$", CITY_KEY, interval))
  water_path <- find_one(DW_DIR, sprintf("^RO_%s_DW_JJA_MULTIYEAR_%s_.*_water\\.tif$", CITY_KEY, interval))
  
  lst_path <- find_one(LST_DIR, sprintf("(?i)^RO_%s_.*LST_HARMONIZED.*JJA.*MULTIYEAR.*%s.*\\.tif$", CITY_KEY, interval))
  
  entropy_paths <- setNames(
    lapply(WINDOWS_M, function(wm) {
      pat1 <- sprintf("(?i)^RO_%s_ENTROPY_%s_W%dm_ALIGNED_TO_SUHI\\.tif$", CITY_KEY, interval, wm)
      f1 <- list.files(ALIGNED_DIR, pattern = pat1, full.names = TRUE)
      if (length(f1) > 0) return(f1[1])
      
      pat2 <- sprintf("(?i)%s.*Hnorm.*%dm.*ALIGNED.*\\.tif$", interval, wm)
      f2 <- list.files(ALIGNED_DIR, pattern = pat2, full.names = TRUE)
      if (length(f2) == 0) stop("Missing entropy aligned file for W", wm, "m in ", ALIGNED_DIR)
      f2[1]
    }),
    paste0("H_", WINDOWS_M)
  )
  
  b <- rast(built_path)
  w <- rast(water_path)
  
  if (!compareGeom(b, w, stopOnError = FALSE)) {
    w <- resample(w, b, method = "bilinear")
  }
  
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
  
  ents <- lapply(entropy_paths, rast)
  for (nm in names(ents)) names(ents[[nm]]) <- nm
  
  for (nm in names(ents)) {
    if (!same.crs(ents[[nm]], lst) || !compareGeom(ents[[nm]], lst, stopOnError = FALSE)) {
      ents[[nm]] <- project(ents[[nm]], lst, method = "bilinear")
    }
  }
  
  AICc_from_AIC <- function(aic, k, n) {
    if (!is.finite(aic) || !is.finite(k) || !is.finite(n)) return(NA_real_)
    if (n <= (k + 1)) return(NA_real_)
    aic + (2 * k * (k + 1)) / (n - k - 1)
  }
  
  safe_lm <- function(formula, data) {
    tryCatch(lm(formula, data = data),
             error = function(e) NULL,
             warning = function(w) invokeRestart("muffleWarning"))
  }
  
  safe_AIC <- function(m) {
    if (is.null(m)) return(NA_real_)
    tryCatch(AIC(m), error = function(e) NA_real_)
  }
  
  safe_k <- function(m) {
    if (is.null(m)) return(NA_real_)
    k <- tryCatch(length(coef(m)), error = function(e) NA_real_)
    if (!is.finite(k)) NA_real_ else k
  }
  
  safe_R2 <- function(m) {
    if (is.null(m)) return(NA_real_)
    tryCatch(summary(m)$r.squared, error = function(e) NA_real_)
  }
  
  core_cand <- ifel((b_lst >= TH_CORE) & (w_lst < TH_WATER), 1, NA)
  p_core <- patches(core_cand, directions = DIRECTIONS)
  fr <- as.data.frame(freq(p_core))
  fr <- fr[!is.na(fr$value) & fr$value != 0, , drop = FALSE]
  if (nrow(fr) == 0) stop("No core candidate pixels found (check TH_CORE/TH_WATER).")
  
  core_id <- fr$value[which.max(fr$count)]
  core <- ifel(p_core == core_id, 1, NA)
  names(core) <- "urban_core"
  
  core_cells <- global(core, "sum", na.rm = TRUE)[1, 1]
  cell_area_m2 <- abs(res(b_lst)[1] * res(b_lst)[2])
  core_area_m2 <- core_cells * cell_area_m2
  r_eq_m <- sqrt(core_area_m2 / pi)
  
  core_dem <- mask(dem_lst, core)
  core_dem_mean <- global(core_dem, "mean", na.rm = TRUE)[1, 1]
  if (!is.finite(core_dem_mean)) stop("Core mean elevation is NA. Check DEM coverage/alignment.")
  
  d <- distance(core)
  
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
  
  make_pts_grid <- function(template_r, grid_m = 300, jitter = TRUE, seed = 1) {
    gridr <- rast(ext(template_r), resolution = grid_m, crs = crs(template_r))
    pts0  <- as.points(gridr, values = FALSE)
    if (!jitter) return(pts0)
    
    xy <- crds(pts0, df = TRUE)
    set.seed(seed)
    xy$x <- xy$x + runif(nrow(xy), -grid_m/2, grid_m/2)
    xy$y <- xy$y + runif(nrow(xy), -grid_m/2, grid_m/2)
    
    e <- ext(template_r)
    xy$x <- pmin(pmax(xy$x, e[1]), e[2])
    xy$y <- pmin(pmax(xy$y, e[3]), e[4])
    
    vect(xy, geom = c("x", "y"), crs = crs(template_r))
  }
  
  PTS_ALL <- make_pts_grid(lst, grid_m = GRID_M, jitter = JITTER_WITHIN_CELL, seed = JITTER_SEED_FIXED)
  
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
  
  make_urban_rural <- function(TH_URBAN, TH_RURAL, BUFFER_M, K_RURAL_MAX, ELEV_TOL_M) {
    RURAL_MAX_M <- as.numeric(K_RURAL_MAX) * r_eq_m
    
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
    names(urban_final) <- "urban"
    
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
    
    list(urban = urban_final, rural = rural, rural_max_m = RURAL_MAX_M)
  }
  
  results_long <- list()
  
  for (i in seq_len(nrow(SENS_GRID))) {
    par <- SENS_GRID[i, ]
    
    urban_rural <- make_urban_rural(
      TH_URBAN    = par$TH_URBAN,
      TH_RURAL    = par$TH_RURAL,
      BUFFER_M    = par$BUFFER_M,
      K_RURAL_MAX = par$K_RURAL_MAX,
      ELEV_TOL_M  = par$ELEV_TOL_M
    )
    
    urban <- urban_rural$urban
    rural <- urban_rural$rural
    
    n_urban_px <- count_cells(urban)
    n_rural_px <- count_cells(rural)
    
    scenario_id <- sprintf(
      "U%.2f_R%.2f_buf%d_k%.1f_tol%d",
      par$TH_URBAN, par$TH_RURAL, par$BUFFER_M, par$K_RURAL_MAX, par$ELEV_TOL_M
    )
    
    if (n_urban_px < 50 || n_rural_px < 50) next
    
    lst_urban <- mask(lst, urban)
    lst_rural <- mask(lst, rural)
    
    med_urban <- g_median(lst_urban)
    med_rural <- g_median(lst_rural)
    
    if (!is.finite(med_urban) || !is.finite(med_rural)) next
    
    suhi_delta_median <- med_urban - med_rural
    
    u_at_pts <- extract(urban, PTS_ALL)[, 2]
    idx_u <- which(is.finite(u_at_pts) & u_at_pts > 0)
    
    if (length(idx_u) < 30) next
    
    pts_u <- PTS_ALL[idx_u]
    
    lst_u <- extract(lst, pts_u)[, 2]
    suhi_u <- lst_u - med_rural
    
    for (nm in names(ents)) {
      h_u <- extract(ents[[nm]], pts_u)[, 2]
      
      ok <- is.finite(suhi_u) & is.finite(h_u)
      n_pts <- sum(ok)
      
      if (n_pts < 30) {
        rho <- NA_real_
        r   <- NA_real_
      } else {
        rho <- suppressWarnings(cor(suhi_u[ok], h_u[ok], method = "spearman"))
        r   <- suppressWarnings(cor(suhi_u[ok], h_u[ok], method = "pearson"))
      }
      
      aic_null <- aic_lin <- aic_quad <- aic_ns3 <- NA_real_
      aicc_null <- aicc_lin <- aicc_quad <- aicc_ns3 <- NA_real_
      best_model <- NA_character_
      best_aicc  <- NA_real_
      delta_aicc_vs_null <- NA_real_
      best_r2 <- NA_real_
      
      if (n_pts >= 30) {
        dfm <- data.frame(y = suhi_u[ok], x = h_u[ok])
        
        m0 <- safe_lm(y ~ 1, data = dfm)
        m1 <- safe_lm(y ~ x, data = dfm)
        m2 <- safe_lm(y ~ poly(x, 2, raw = TRUE), data = dfm)
        m3 <- safe_lm(y ~ splines::ns(x, df = 3), data = dfm)
        
        aic_null <- safe_AIC(m0); k0 <- safe_k(m0); aicc_null <- AICc_from_AIC(aic_null, k0, n_pts)
        aic_lin  <- safe_AIC(m1); k1 <- safe_k(m1); aicc_lin  <- AICc_from_AIC(aic_lin,  k1, n_pts)
        aic_quad <- safe_AIC(m2); k2 <- safe_k(m2); aicc_quad <- AICc_from_AIC(aic_quad, k2, n_pts)
        aic_ns3  <- safe_AIC(m3); k3 <- safe_k(m3); aicc_ns3  <- AICc_from_AIC(aic_ns3,  k3, n_pts)
        
        aicc_vec <- c(null = aicc_null, lin = aicc_lin, quad = aicc_quad, ns3 = aicc_ns3)
        aicc_vec <- aicc_vec[is.finite(aicc_vec)]
        
        if (length(aicc_vec) > 0) {
          best_model <- names(which.min(aicc_vec))
          best_aicc  <- min(aicc_vec)
          
          if (is.finite(aicc_null)) delta_aicc_vs_null <- best_aicc - aicc_null
          
          best_fit <- switch(best_model,
                             null = m0,
                             lin  = m1,
                             quad = m2,
                             ns3  = m3)
          if (best_model == "null") {
            best_r2 <- 0
          } else {
            best_r2 <- safe_R2(best_fit)
          }
        }
      }
      
      results_long[[length(results_long) + 1]] <- data.frame(
        city = CITY_KEY,
        interval = interval,
        scenario_id = scenario_id,
        
        TH_URBAN = par$TH_URBAN,
        TH_RURAL = par$TH_RURAL,
        BUFFER_M = par$BUFFER_M,
        K_RURAL_MAX = par$K_RURAL_MAX,
        RURAL_MAX_M = urban_rural$rural_max_m,
        ELEV_TOL_M = par$ELEV_TOL_M,
        
        n_urban_px = n_urban_px,
        n_rural_px = n_rural_px,
        
        lst_urban_median = med_urban,
        lst_rural_median = med_rural,
        suhi_delta_median_C = suhi_delta_median,
        
        window = nm,
        n_points = n_pts,
        spearman_rho = rho,
        pearson_r = r,
        
        aic_null = aic_null,
        aic_lin  = aic_lin,
        aic_quad = aic_quad,
        aic_ns3  = aic_ns3,
        
        aicc_null = aicc_null,
        aicc_lin  = aicc_lin,
        aicc_quad = aicc_quad,
        aicc_ns3  = aicc_ns3,
        
        best_model = best_model,
        best_aicc  = best_aicc,
        delta_aicc_vs_null = delta_aicc_vs_null,
        best_r2 = best_r2,
        
        row.names = NULL
      )
    }
    
    rm(lst_urban, lst_rural)
    gc()
  }
  
  res_long <- if (length(results_long)) do.call(rbind, results_long) else data.frame()
  
  out_long <- file.path(OUT_DIR, sprintf("RO_%s_ sensitivity_results_LONG_%s.csv", CITY_KEY, interval))
  write.csv(res_long, out_long, row.names = FALSE)
}