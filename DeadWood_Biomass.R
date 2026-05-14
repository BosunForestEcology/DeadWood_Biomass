# spatialExtent field omitted: removed from SpaDES.core API in version >= 3.0
defineModule(sim, list(
  name        = "DeadWood_Biomass",
  description = "Translates snag and DWD decay class inventories into pixel-level biomass
                 estimates (Mg ha-1) using species- and pool-specific density reduction
                 factors (DRF) from Paper 2 Appendix D.",
  keywords    = c("dead wood", "biomass", "density reduction factor", "carbon"),
  authors     = structure(list(list(given = "Thomson", family = "Harris",
                                    role = c("aut", "cre"),
                                    email = "BosunForestEcology@gmail.com", comment = NULL)),
                           class = "person"),
  childModules = character(0),
  version     = list(DeadWood_Biomass = "0.0.1"),
  timeframe   = as.POSIXlt(c(NA, NA)),
  timeunit    = "year",
  citation    = list(),
  documentation = list(),
  reqdPkgs    = list("data.table", "terra", "ggplot2", "SpaDES.core (>= 3.0.0)"),
  parameters  = bindrows(
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    desc = "Simulation time for first plot. NA = no plots."),
    defineParameter(".plotInterval",    "numeric",  5, NA, NA,
                    desc = "Interval between plots (years). Defaults to 5 to match transition timestep.")
  ),
  inputObjects = bindrows(
    expectsInput("snagTable", "data.table",
                 desc = "Current snag inventory from snagDecay."),
    expectsInput("DWDTable", "data.table",
                 desc = "Current DWD inventory from DWDDecay."),
    expectsInput("studyAreaRaster", "SpatRaster",
                 desc = "Template raster defining pixel grid, CRS, and resolution.")
  ),
  outputObjects = bindrows(
    createsOutput("snagBiomass_Mg_ha", "SpatRaster",
                  desc = "Pixel-level snag biomass at death (Mg ha-1), summed across decay classes. No DRF applied."),
    createsOutput("DWDBiomass_Mg_ha", "SpatRaster",
                  desc = "Pixel-level DWD biomass at death (Mg ha-1), summed across decay classes. No DRF applied."),
    createsOutput("snagHistory", "SpatRaster",
                  desc = "Multi-layer raster of snag biomass snapshots at each plot interval (all species combined)."),
    createsOutput("DWDHistory", "SpatRaster",
                  desc = "Multi-layer raster of DWD biomass snapshots at each plot interval (all species combined)."),
    createsOutput("snagHistoryBySpecies", "list",
                  desc = "Named list of multi-layer SpatRasters of snag biomass snapshots per species."),
    createsOutput("DWDHistoryBySpecies", "list",
                  desc = "Named list of multi-layer SpatRasters of DWD biomass snapshots per species.")
  )
))

doEvent.DeadWood_Biomass <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, start(sim) + 5, "DeadWood_Biomass", "transition", eventPriority = 4)
      if (!is.na(P(sim)$.plotInitialTime))
        sim <- scheduleEvent(sim, P(sim)$.plotInitialTime,
                             "DeadWood_Biomass", "plot", eventPriority = 5)
    },
    transition = {
      sim <- Transition(sim)
      sim <- scheduleEvent(sim, time(sim) + 5, "DeadWood_Biomass", "transition", eventPriority = 4)
    },
    plot = {
      sim <- Snapshot(sim)
      nextPlot <- time(sim) + P(sim)$.plotInterval
      if (nextPlot <= end(sim))
        sim <- scheduleEvent(sim, nextPlot, "DeadWood_Biomass", "plot", eventPriority = 5)
    },
    warning(paste("Undefined event type:", eventType, "in module deadWoodBiomass"))
  )
  return(invisible(sim))
}

Init <- function(sim) {
  # terra::values<- deep-copies before writing, so studyAreaRaster is not modified
  sim$snagBiomass_Mg_ha <- sim$studyAreaRaster
  terra::values(sim$snagBiomass_Mg_ha) <- NA_real_
  sim$DWDBiomass_Mg_ha  <- sim$studyAreaRaster
  terra::values(sim$DWDBiomass_Mg_ha)  <- NA_real_
  sim$snagHistory          <- NULL
  sim$DWDHistory           <- NULL
  sim$snagHistoryBySpecies <- list()
  sim$DWDHistoryBySpecies  <- list()
  return(invisible(sim))
}

Transition <- function(sim) {
  # Sum initBiomass (biomass at death) by pixel — no DRF applied
  snagByPixel <- sim$snagTable[, .(value = sum(initBiomass, na.rm = TRUE)), by = pixelID]
  sim$snagBiomass_Mg_ha <- pixelValuesToRaster(snagByPixel, sim$studyAreaRaster)

  DWDbyPixel <- sim$DWDTable[, .(value = sum(initBiomass, na.rm = TRUE)), by = pixelID]
  sim$DWDBiomass_Mg_ha <- pixelValuesToRaster(DWDbyPixel, sim$studyAreaRaster)

  return(invisible(sim))
}

Snapshot <- function(sim) {
  yr        <- as.integer(time(sim))
  snagLayer <- sim$snagBiomass_Mg_ha
  DWDLayer  <- sim$DWDBiomass_Mg_ha
  names(snagLayer) <- paste0("yr", yr)
  names(DWDLayer)  <- paste0("yr", yr)
  sim$snagHistory <- if (is.null(sim$snagHistory)) snagLayer else c(sim$snagHistory, snagLayer)
  sim$DWDHistory  <- if (is.null(sim$DWDHistory))  DWDLayer  else c(sim$DWDHistory,  DWDLayer)

  allSp <- unique(c(sim$snagTable$species, sim$DWDTable$species))
  for (sp in allSp) {
    snagSp      <- sim$snagTable[species == sp, .(value = sum(initBiomass, na.rm = TRUE)), by = pixelID]
    spSnagLayer <- pixelValuesToRaster(snagSp, sim$studyAreaRaster)
    names(spSnagLayer) <- paste0("yr", yr)
    sim$snagHistoryBySpecies[[sp]] <-
      if (is.null(sim$snagHistoryBySpecies[[sp]])) spSnagLayer else c(sim$snagHistoryBySpecies[[sp]], spSnagLayer)

    DWDsp      <- sim$DWDTable[species == sp, .(value = sum(initBiomass, na.rm = TRUE)), by = pixelID]
    spDWDLayer <- pixelValuesToRaster(DWDsp, sim$studyAreaRaster)
    names(spDWDLayer) <- paste0("yr", yr)
    sim$DWDHistoryBySpecies[[sp]] <-
      if (is.null(sim$DWDHistoryBySpecies[[sp]])) spDWDLayer else c(sim$DWDHistoryBySpecies[[sp]], spDWDLayer)
  }

  return(invisible(sim))
}
