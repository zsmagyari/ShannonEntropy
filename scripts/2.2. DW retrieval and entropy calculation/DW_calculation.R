library(terra)

CITY_KEY <- "CLUJ_NAPOCA"
YEARS <- "2021-2025"

in_dir  <- file.path("DW",CITY_KEY)         
out_dir <- file.path("DW",CITY_KEY,"entropy_no_qc") 
tmp_dir <- "F:/RTMP"         

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

WINDOW_M <- c(90, 150, 300, 600)
LOG_BASE <- "e"  
EXPORT_NORMALIZED <- TRUE

terraOptions(
  tempdir = tmp_dir,
  todisk  = TRUE,
  memfrac = 0.60,
  progress = 1
)

WOPT_FLOAT <- list(
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6", "TILED=YES", "BIGTIFF=YES"),
  datatype = "FLT4S"
)

bands_canon <- c(
  "water", "trees", "grass", "flooded_vegetation", "crops",
  "shrub_and_scrub", "built", "bare", "snow_and_ice"
)

band_pattern <- list(
  water = "water",
  trees = "(trees|tree)",
  grass = "(grass|frass)",
  flooded_vegetation = "(flooded_vegetation|fooded_vegetation|flooded-vegetation)",
  crops = "crops",
  shrub_and_scrub = "shrub_and_scrub",
  built = "built",
  bare = "bare",
  snow_and_ice = "snow_and_ice"
)

log_fun <- function(x) {
  if (LOG_BASE == "2") return(log(x, base = 2))
  log(x)
}

entropy_vec <- function(v) {
  if (all(is.na(v))) return(NA_real_)
  v[is.na(v)] <- 0
  v[v < 0] <- 0
  s <- sum(v)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  v <- v / s
  v <- v[v > 0]
  -sum(v * log_fun(v))
}

window_ncells <- function(win_m, res_m) {
  n <- as.integer(round(win_m / res_m))
  if (n < 1) n <- 1L
  if (n %% 2 == 0) n <- n + 1L
  n
}

find_band_file <- function(interval, band_key) {
  pat_band <- band_pattern[[band_key]]
  pat <- sprintf("^RO_%s_DW_JJA_MULTIYEAR_%s_.*_%s\\.tif$", CITY_KEY, interval, pat_band)
  f <- list.files(in_dir, pattern = pat, full.names = TRUE, ignore.case = TRUE)
  if (length(f) == 0) stop(sprintf("Missing file: band '%s'", band_key))
  f[[1]]
}

safe_tmp <- function(prefix, interval, win_m = NA_integer_, ext = ".tif") {
  if (is.na(win_m)) {
    file.path(tmp_dir, sprintf("TMP_%s_%s_%s%s", CITY_KEY, prefix, interval, ext))
  } else {
    file.path(tmp_dir, sprintf("TMP_%s_%s_%s_win%dm%s", CITY_KEY, prefix, interval, win_m, ext))
  }
}

calc_entropy_term <- function(x) {
  base_val <- if(LOG_BASE == "2") 2 else exp(1)
  lx <- log(x, base = base_val)
  res <- x * lx
  res[!is.finite(res)] <- 0
  return(res)
}

files <- vapply(bands_canon, function(b) find_band_file(YEARS, b), character(1))
r_list <- lapply(files, rast)
template <- r_list[[1]]

for (i in seq_along(r_list)) {
  if (!compareGeom(template, r_list[[i]], stopOnError = FALSE)) {
    r_list[[i]] <- resample(r_list[[i]], template, method = "bilinear")
  }
}

P <- rast(r_list)
names(P) <- bands_canon

sumP_fn <- safe_tmp("sumP", YEARS)
sumP <- app(P, fun = function(v) sum(v, na.rm = TRUE),
            filename = sumP_fn, overwrite = TRUE, wopt = WOPT_FLOAT)

Pn_fn <- safe_tmp("Pn", YEARS)
Pn <- (P / sumP)
Pn <- mask(Pn, sumP > 0)
Pn <- writeRaster(Pn, filename = Pn_fn, overwrite = TRUE, wopt = WOPT_FLOAT)

H0_fn <- safe_tmp("H0_pixel", YEARS)
H0 <- app(Pn, fun = entropy_vec,
          filename = H0_fn, overwrite = TRUE, wopt = WOPT_FLOAT)

if (EXPORT_NORMALIZED) {
  H0n <- H0 / log_fun(length(bands_canon))
  out0 <- file.path(out_dir, sprintf("RO_%s_DW_Hnorm_JJA_%s_pixel.tif", CITY_KEY, YEARS))
  writeRaster(H0n, out0, overwrite = TRUE, wopt = WOPT_FLOAT)
  rm(H0n)
} else {
  out0 <- file.path(out_dir, sprintf("RO_%s_DW_H_JJA_%s_pixel.tif", CITY_KEY, YEARS))
  writeRaster(H0, out0, overwrite = TRUE, wopt = WOPT_FLOAT)
}
rm(H0)

res_m <- res(Pn)[1]

for (win_m in WINDOW_M) {
  n <- window_ncells(win_m, res_m)
  w <- matrix(1, nrow = n, ncol = n)
  message(sprintf("\nProcessing Window: %d m...", win_m))
  
  H_accum_fn <- safe_tmp("H_accum_init", YEARS, win_m)
  H_accum <- init(Pn[[1]], 0, filename=H_accum_fn, overwrite=TRUE, wopt=WOPT_FLOAT)
  
  for (b_idx in 1:nlyr(Pn)) {
    single_band <- Pn[[b_idx]]
    
    band_mean_fn <- safe_tmp(paste0("band", b_idx, "_mean"), YEARS, win_m)
    band_mean <- focal(single_band, w = w, fun = "mean", na.rm = TRUE, 
                       filename = band_mean_fn, overwrite = TRUE, wopt = WOPT_FLOAT)
    
    term <- calc_entropy_term(band_mean)
    
    H_new_fn <- safe_tmp(paste0("H_accum_step", b_idx), YEARS, win_m)
    H_next <- H_accum - term
    H_next <- writeRaster(H_next, filename = H_new_fn, overwrite = TRUE, wopt = WOPT_FLOAT)
    
    prev_fn <- sources(H_accum)
    rm(H_accum, band_mean, term, single_band)
    gc()
    if (file.exists(prev_fn)) file.remove(prev_fn)
    if (file.exists(band_mean_fn)) file.remove(band_mean_fn)
    
    H_accum <- H_next
  }
  
  if (EXPORT_NORMALIZED) {
    base_val <- if(LOG_BASE == "2") 2 else exp(1)
    max_ent <- log(nlyr(Pn), base = base_val)
    Hw <- H_accum / max_ent
    Hw <- clamp(Hw, lower=0, upper=1)
    
    outH <- file.path(out_dir, sprintf("RO_%s_DW_Hnorm_JJA_%s_win%dm.tif", CITY_KEY, YEARS, win_m))
    writeRaster(Hw, outH, overwrite = TRUE, wopt = WOPT_FLOAT)
  } else {
    outH <- file.path(out_dir, sprintf("RO_%s_DW_H_JJA_%s_win%dm.tif", CITY_KEY, YEARS, win_m))
    writeRaster(H_accum, outH, overwrite = TRUE, wopt = WOPT_FLOAT)
  }
  
  rm(H_accum)
  gc()
  tmp_pattern <- sprintf("TMP_%s_.*_win%dm.*", CITY_KEY, win_m)
  unlink(list.files(tmp_dir, pattern = tmp_pattern, full.names = TRUE))
}
