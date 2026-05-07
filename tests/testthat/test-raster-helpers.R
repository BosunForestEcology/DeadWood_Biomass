library(terra)
library(data.table)
library(testthat)

# This test file lives at: modules/deadWoodBiomass/tests/testthat/
# The project root (containing modules/) is 4 levels up.
.projRoot <- normalizePath(file.path(getwd(), "..", "..", "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(.projRoot, "modules"))) {
  .projRoot <- getwd()
}
source(file.path(.projRoot, "modules", "deadWoodBiomass", "R", "raster-helpers.R"))

test_that("pixelValuesToRaster assigns values to correct pixels", {
  template <- terra::rast(nrows = 3, ncols = 3,
                          xmin = 0, xmax = 3, ymin = 0, ymax = 3,
                          crs = "EPSG:4326")
  dt <- data.table(pixelID = c(1L, 5L, 9L), value = c(10.0, 20.0, 30.0))
  r <- pixelValuesToRaster(dt, template)
  expect_equal(unname(terra::values(r)[1, 1]),  10.0)
  expect_equal(unname(terra::values(r)[5, 1]),  20.0)
  expect_equal(unname(terra::values(r)[9, 1]),  30.0)
  expect_true(is.na(terra::values(r)[2, 1]))
})

test_that("pixelValuesToRaster returns NA for unspecified pixels", {
  template <- terra::rast(nrows = 2, ncols = 2,
                          xmin = 0, xmax = 2, ymin = 0, ymax = 2,
                          crs = "EPSG:4326")
  dt <- data.table(pixelID = 1L, value = 99.0)
  r <- pixelValuesToRaster(dt, template)
  expect_equal(sum(is.na(terra::values(r))), 3L)
})

test_that("pixelValuesToRaster handles empty input data.table", {
  template <- terra::rast(nrows = 2, ncols = 2,
                          xmin = 0, xmax = 2, ymin = 0, ymax = 2,
                          crs = "EPSG:4326")
  dt <- data.table(pixelID = integer(), value = numeric())
  r <- pixelValuesToRaster(dt, template)
  expect_true(all(is.na(terra::values(r))))
})
