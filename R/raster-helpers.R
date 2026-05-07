# Map pixel IDs and values from a data.table onto a copy of templateRaster.
# dt must have columns pixelID (integer) and value (numeric).
# Returns a SpatRaster with NA for all pixels not in dt.
pixelValuesToRaster <- function(dt, templateRaster) {
  r <- templateRaster
  terra::values(r) <- NA_real_
  if (nrow(dt) > 0) {
    r[dt$pixelID] <- dt$value
  }
  r
}
