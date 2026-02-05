library(terra)

rasters_root <- "D:/11_COHERENCE/1_DATA_SERIES/UNZIPPED"

corr_files <- list.files(
  rasters_root,
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE,
  pattern = "\\.(tif|tiff)$"
)

corr_files <- corr_files[grepl("_corr", basename(corr_files), ignore.case = TRUE)]
if (length(corr_files) == 0) stop("No *_corr TIFF files found under: ", rasters_root)

# ---- extract date between first "_" and "T" ----
get_date <- function(x) {
  name <- basename(x)
  sub("^.*?_([0-9]{8})T.*$", "\\1", name)
}

dates_chr <- get_date(corr_files)

# Convert to Date (safer than character sorting)
dates <- as.Date(dates_chr, format = "%Y%m%d")

# Order files by date
corr_files <- corr_files[order(dates)]


# Read
r_list <- lapply(corr_files, rast)
ref <- r_list[[1]]

# Align to ref (CRS + grid)
r_list2 <- lapply(r_list, function(r) {
  if (!same.crs(r, ref)) r <- project(r, ref)
  resample(r, ref, method = "bilinear")
})

# Crop all to common overlap
common_ext <- Reduce(intersect, lapply(r_list2, ext))
r_list3 <- lapply(r_list2, function(r) crop(r, common_ext))

# Stack and write
corr_stack <- rast(r_list3)

writeRaster(
  corr_stack *  10000, # scale factor
  filename = file.path(rasters_root, "corr_stack"),
  filetype = "ENVI",  
  datatype = "INT2S",
  gdal = c("COMPRESS=DEFLATE"),
  overwrite = TRUE
)