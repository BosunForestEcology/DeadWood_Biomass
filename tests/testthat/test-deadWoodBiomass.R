library(SpaDES.core)
library(data.table)
library(terra)
library(testthat)

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
    params  = list(DeadWood_Biomass = list()),
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

test_that("deadWoodBiomass transition computes snag biomass as initBiomass sum", {
  snagTable <- data.table(
    pixelID = 1L, species = "Pinus strobus", DC = 1L, ageInDC = 0L, initBiomass = 12.0
  )
  sim <- testInit(
    "DeadWood_Biomass",
    times   = list(start = 0, end = 5),
    params  = list(DeadWood_Biomass = list()),
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

test_that("deadWoodBiomass transition computes DWD biomass as initBiomass sum", {
  DWDTable <- data.table(
    pixelID = 2L, species = "Pinus strobus", DC = 3L, ageInDC = 1L, initBiomass = 10.0
  )
  sim <- testInit(
    "DeadWood_Biomass",
    times   = list(start = 0, end = 5),
    params  = list(DeadWood_Biomass = list()),
    objects = list(
      snagTable       = data.table::copy(emptySnagTable),
      DWDTable        = DWDTable,
      studyAreaRaster = templateRaster
    )
  )
  sim <- spades(sim, events = c("init", "transition"))
  expect_equal(unname(terra::values(sim$DWDBiomass_Mg_ha)[2, 1]), 10.0)
})

test_that("deadWoodBiomass transition aggregates multiple records per pixel", {
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
    params  = list(DeadWood_Biomass = list()),
    objects = list(
      snagTable       = snagTable,
      DWDTable        = data.table::copy(emptySnagTable),
      studyAreaRaster = templateRaster
    )
  )
  sim <- spades(sim, events = c("init", "transition"))
  expect_equal(unname(terra::values(sim$snagBiomass_Mg_ha)[3, 1]), 13.0)
})
