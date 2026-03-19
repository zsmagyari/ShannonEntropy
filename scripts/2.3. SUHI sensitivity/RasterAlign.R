library(terra)

CITY_KEY <- "TIMISOARA"
interval    <- "2021-2025"

WINDOWS_M <- c(90, 150, 300, 600)

suhi_dir    <- file.path("SUHI",CITY_KEY)
entropy_dir <- file.path("DW",CITY_KEY,"entropy")
out_dir     <- file.path("ANALYSIS","aligned",CITY_KEY)

tmp_dir <- "F:/RTMP"              
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

terraOptions(tempdir = tmp_dir, todisk = TRUE, memfrac = 0.6, progress = 1)

WOPT <- list(
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6", "TILED=YES", "BIGTIFF=YES")
)

find_suhi_file <- function(interval) {
  pat <- sprintf("(?i)RO_%s_.*SUHI.*%s.*\\.tif$", CITY_KEY, interval)
  f <- list.files(suhi_dir, pattern = pat, full.names = TRUE)
  if (length(f) == 0) stop("Missing SUHI file: year=", interval, " (pattern: ", pat, ")")
  if (length(f) > 1) message("Many SUHI files (", interval, "), first will be used: ", basename(f[1]))
  f[1]
}

find_entropy_file <- function(interval, win_m) {
  pat <- sprintf("(?i)RO_%s_.*Hnorm.*%s.*%dm.*\\.tif$", CITY_KEY, interval, win_m)
  f <- list.files(entropy_dir, pattern = pat, full.names = TRUE)
  if (length(f) == 0) stop("Missing entropy file: year=", interval, " win=", win_m, "m (pattern: ", pat, ")")
  if (length(f) > 1) message("Many entropy files  (", interval, ", ", win_m, "m), first will be used: ", basename(f[1]))
  f[1]
}

  suhi_path <- find_suhi_file(interval)
  suhi <- rast(suhi_path)
  
  suhi_out <- file.path(out_dir, sprintf("RO_%s_SUHI_%s_ALIGNED_TEMPLATE.tif", CITY_KEY, interval))
  writeRaster(suhi, suhi_out, overwrite = TRUE, wopt = WOPT, datatype = "FLT4S")
  
  message("SUHI template: ", basename(suhi_path))
  message("  CRS:  ", crs(suhi))
  message("  Res:  ", paste(res(suhi), collapse = " x "))
  message("  Ext:  ", paste(as.vector(ext(suhi)), collapse = ", "))
  
  for (win in WINDOWS_M) {
    message("  -> Entropy window: ", win, " m")
    
    ent_path <- find_entropy_file(interval, win)
    ent <- rast(ent_path)
    
    if (!same.crs(ent, suhi)) {
      ent_al <- project(ent, suhi, method = "bilinear")
    } else if (!compareGeom(ent, suhi, stopOnError = FALSE)) {
      ent_al <- resample(ent, suhi, method = "bilinear")
    } else {
      ent_al <- ent
    }

    out_ent <- file.path(out_dir, sprintf("RO_%s_ENTROPY_%s_W%dm_ALIGNED_TO_SUHI.tif", CITY_KEY, interval, win))
    writeRaster(ent_al, out_ent, overwrite = TRUE, wopt = WOPT, datatype = "FLT4S")    
    rm(ent, ent_al); gc()
  }
  
  rm(suhi); gc()
