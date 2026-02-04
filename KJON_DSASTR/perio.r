#### Script to calculate periodogram of time series extracted from raster data at given point locations
library(terra)
library(ncdf4)

#' Read multiband raster data
#'
#' @param input_path Path to input raster file (GeoTiff or NetCDF)
#' @return SpatRaster object with multiple layers (bands)
read_raster <- function(input_path) {
  # Check file extension to handle NetCDF specially

if (grepl("\\.nc$", input_path, ignore.case = TRUE)) {
    raster <- rast(input_path)
  } else {
    raster <- rast(input_path)
  }
  return(raster)
}

#' Calculate periodogram for a single time series
#'
#' @param ts Numeric vector representing a time series
#' @return List with frequency and psd (power spectral density)
calculate_single_periodogram <- function(ts) {
  n <- length(ts)
  # Remove mean (detrend)
  ts_centered <- ts - mean(ts, na.rm = TRUE)
  
  # FFT
  fft_result <- fft(ts_centered)
  
  # Power spectral density (one-sided)
  psd <- (Mod(fft_result)^2) / n
  
  # Frequencies (normalized, 0 to 0.5)
  freq <- seq(0, 0.5, length.out = floor(n / 2) + 1)
  
  # Take only positive frequencies (one-sided spectrum)
  n_freq <- length(freq)
  psd <- psd[1:n_freq]
  
  # Scale for one-sided (multiply by 2, except DC and Nyquist)
  psd[2:(n_freq - 1)] <- 2 * psd[2:(n_freq - 1)]
  
  return(list(freq = freq, psd = psd))
}

#' Calculate periodogram for each pixel's time series
#'
#' @param raster SpatRaster with multiple layers (time series per pixel)
#' @return List with psd_raster (SpatRaster with PSD values) and freq (frequency bins)
calculate_periodogram_raster <- function(raster) {
  n_layers <- nlyr(raster)
  n_rows <- nrow(raster)
  n_cols <- ncol(raster)
  
  # Get frequency bins from sample calculation
  sample_ts <- rep(0, n_layers)
  sample_result <- calculate_single_periodogram(sample_ts)
  freq <- sample_result$freq
  n_freq <- length(freq)
  
  message(sprintf("Processing %d x %d pixels with %d time steps -> %d frequencies",
                  n_rows, n_cols, n_layers, n_freq))
  
  # Extract all values as matrix (pixels x layers)
  vals <- values(raster)  # matrix: n_cells x n_layers
  n_cells <- nrow(vals)
  
  # Prepare output matrix
  psd_matrix <- matrix(NA_real_, nrow = n_cells, ncol = n_freq)
  
  # Calculate periodogram for each pixel
  pb <- txtProgressBar(min = 0, max = n_cells, style = 3)
  for (i in seq_len(n_cells)) {
    ts <- vals[i, ]
    if (!any(is.na(ts))) {
      result <- calculate_single_periodogram(ts)
      psd_matrix[i, ] <- result$psd
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Create output raster with PSD values (one layer per frequency)
  psd_raster <- rast(raster, nlyrs = n_freq)
  values(psd_raster) <- psd_matrix
  
  # Set layer names
  names(psd_raster) <- paste0("freq_", round(freq, 4))
  
  # Copy CRS from original
  crs(psd_raster) <- crs(raster)
  ext(psd_raster) <- ext(raster)
  
  return(list(psd_raster = psd_raster, freq = freq))
}

#' Save PSD raster to GeoTiff file
#'
#' @param output_path Path to output raster file
#' @param psd_raster SpatRaster with PSD values
save_raster <- function(output_path, psd_raster) {
  writeRaster(psd_raster, output_path,
              overwrite = TRUE,
              datatype = "FLT4S",
              gdal = c("COMPRESS=LZW"))
}

#' Main function
main <- function() {
  input_path <- "E:/SAR_UVa/coherence_stack.nc"
  output_path <- "E:/SAR_UVa/output_psd_raster.tif"
  periods_path <- "E:/SAR_UVa/output_periods.txt"
  
  # Read raster data
  message("Reading raster data...")
  data <- read_raster(input_path)
  
  # Calculate periodogram
  message("Calculating periodogram...")
  result <- calculate_periodogram_raster(data)
  psd_raster <- result$psd_raster
  freq <- result$freq
  
  # Save PSD raster
  message("Saving PSD raster...")
  save_raster(output_path, psd_raster)
  
  # Save periods (1/frequency)
  message("Saving periods...")
  # Avoid division by zero for DC component
  periods <- ifelse(freq == 0, Inf, 1 / freq)
  write.table(periods, periods_path,
              row.names = FALSE, col.names = FALSE)
  
  message("Done!")
}

# Run main
if (!interactive()) {
  main()
}
