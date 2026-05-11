library(SpaDES.core)
library(data.table)
library(terra)
library(testthat)

DRFLookup_test <- data.table(
  species = "Pinus strobus",
  pool    = rep(c("snag", "DWD"), each = 5),
  DC      = rep(1:5, times = 2),
  DRF     = c(1.000, 0.841, 0.706, 0.543, 0.382,
              1.000, 0.783, 0.614, 0.418, 0.251)
)

templateRaster <- terra::rast(nrows = 3, ncols = 3,
                               xmin = 0, xmax = 3, ymin = 0, ymax = 3,
                               crs = "EPSG:4326")

emptySnagTable <- data.table(
  pixelID = integer(), species = character(),
  DC = integer(), ageInDC = integer(), initBiomass = numeric()
)

# testInit is not available in SpaDES.core >= 3.x; use simInit directly.
# This test file lives at: modules/DeadWood_Biomass/tests/testthat/
# The project root (containing modules/) is 4 levels up.
.projRoot <- normalizePath(file.path(getwd(), "..", "..", "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(.projRoot, "modules"))) {
  .projRoot <- getwd()
}

testInit <- function(moduleName, params, objects, times = list(start = 0, end = 10)) {
  simInit(
    times   = times,
    modules = list(moduleName),
    params  = params,
    objects = objects,
    paths   = list(modulePath = file.path(.projRoot, "modules"))
  )
}

test_that("deadWoodBiomass init creates snagBiomass and DWDBiomass rasters", {
  sim <- testInit(
    "DeadWood_Biomass",
    params  = list(DeadWood_Biomass = list(DRFLookup = DRFLookup_test)),
    objects = list(
      snagTable       = data.table::copy(emptySnagTable),
      DWDTable        = data.table::copy(emptySnagTable),
      studyAreaRaster = templateRaster
    )
  )
  sim <- spades(sim, events = "init")
  expect_s4_class(sim$snagBiomass_Mg_ha, "SpatRaster")
  expect_s4_class(sim$DWDBiomass_Mg_ha,  "SpatRaster")
  expect_true(all(is.na(terra::values(sim$snagBiomass_Mg_ha))))
  expect_true(all(is.na(terra::values(sim$DWDBiomass_Mg_ha))))
})

test_that("deadWoodBiomass transition computes snag biomass correctly for DC1", {
  # DC1 snag: DRF = 1.0, so currentBiomass == initBiomass
  snagTable <- data.table(
    pixelID = 1L, species = "Pinus strobus", DC = 1L, ageInDC = 0L, initBiomass = 12.0
  )
  sim <- testInit(
    "DeadWood_Biomass",
    times   = list(start = 0, end = 5),
    params  = list(DeadWood_Biomass = list(DRFLookup = DRFLookup_test)),
    objects = list(
      snagTable       = snagTable,
      DWDTable        = data.table::copy(emptySnagTable),
      studyAreaRaster = templateRaster
    )
  )
  sim <- spades(sim, events = c("init", "transition"))
  expect_equal(unname(terra::values(sim$snagBiomass_Mg_ha)[1, 1]), 12.0)
  expect_true(all(is.na(terra::values(sim$DWDBiomass_Mg_ha))))
})

test_that("deadWoodBiomass transition applies DRF correctly for DC3 DWD", {
  # DC3 DWD DRF = 0.614; initBiomass = 10 -> currentBiomass = 6.14
  DWDTable <- data.table(
    pixelID = 2L, species = "Pinus strobus", DC = 3L, ageInDC = 1L, initBiomass = 10.0
  )
  sim <- testInit(
    "DeadWood_Biomass",
    times   = list(start = 0, end = 5),
    params  = list(DeadWood_Biomass = list(DRFLookup = DRFLookup_test)),
    objects = list(
      snagTable       = data.table::copy(emptySnagTable),
      DWDTable        = DWDTable,
      studyAreaRaster = templateRaster
    )
  )
  sim <- spades(sim, events = c("init", "transition"))
  expect_equal(unname(terra::values(sim$DWDBiomass_Mg_ha)[2, 1]), 10.0 * 0.614, tolerance = 1e-6)
})

test_that("deadWoodBiomass transition aggregates multiple records per pixel", {
  # Two DC2 snags in pixel 3: 5 + 8 = 13 Mg/ha, each x DRF 0.841 = 10.933
  snagTable <- data.table(
    pixelID     = c(3L, 3L),
    species     = "Pinus strobus",
    DC          = c(2L, 2L),
    ageInDC     = c(0L, 0L),
    initBiomass = c(5.0, 8.0)
  )
  sim <- testInit(
    "DeadWood_Biomass",
    times   = list(start = 0, end = 5),
    params  = list(DeadWood_Biomass = list(DRFLookup = DRFLookup_test)),
    objects = list(
      snagTable       = snagTable,
      DWDTable        = data.table::copy(emptySnagTable),
      studyAreaRaster = templateRaster
    )
  )
  sim <- spades(sim, events = c("init", "transition"))
  expect_equal(unname(terra::values(sim$snagBiomass_Mg_ha)[3, 1]), (5.0 + 8.0) * 0.841, tolerance = 1e-6)
})
