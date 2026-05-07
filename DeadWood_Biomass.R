# spatialExtent field omitted: removed from SpaDES.core API in version >= 3.0
defineModule(sim, list(
  name        = "DeadWood_Biomass",
  description = "Translates snag and DWD decay class inventories into pixel-level biomass
                 estimates (Mg ha-1) using species- and pool-specific density reduction
                 factors (DRF) from Paper 2 Appendix D.",
  keywords    = c("dead wood", "biomass", "density reduction factor", "carbon"),
  authors     = structure(list(list(given = "First", family = "Last",
                                    role = c("aut", "cre"),
                                    email = "email@example.com", comment = NULL)),
                           class = "person"),
  childModules = character(0),
  version     = list(DeadWood_Biomass = "0.0.1"),
  timeframe   = as.POSIXlt(c(NA, NA)),
  timeunit    = "year",
  citation    = list(),
  documentation = list(),
  reqdPkgs    = list("data.table", "terra", "ggplot2", "SpaDES.core (>= 3.0.0)"),
  parameters  = bindrows(
    defineParameter("DRFLookup", "data.table",
                    data.table::data.table(species = character(), pool = character(),
                                           DC = integer(), DRF = numeric()),
                    NA, NA,
                    desc = "Density reduction factor lookup: species, pool, DC, DRF."),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    desc = "Simulation time for first plot. NA = no plots."),
    defineParameter(".plotInterval",    "numeric",  1, NA, NA,
                    desc = "Interval between plots.")
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
                  desc = "Pixel-level snag biomass (Mg ha-1)."),
    createsOutput("DWDBiomass_Mg_ha", "SpatRaster",
                  desc = "Pixel-level DWD biomass (Mg ha-1)."),
    createsOutput("snagHistory", "SpatRaster",
                  desc = "Multi-layer raster of snag biomass snapshots at each plot interval."),
    createsOutput("DWDHistory", "SpatRaster",
                  desc = "Multi-layer raster of DWD biomass snapshots at each plot interval.")
  )
))

doEvent.DeadWood_Biomass <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, start(sim) + 1, "DeadWood_Biomass", "annual", eventPriority = 4)
      if (!is.na(P(sim)$.plotInitialTime))
        sim <- scheduleEvent(sim, P(sim)$.plotInitialTime,
                             "DeadWood_Biomass", "plot", eventPriority = 5)
    },
    annual = {
      sim <- Annual(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "DeadWood_Biomass", "annual", eventPriority = 4)
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
  if (nrow(P(sim)$DRFLookup) == 0L)
    stop("DRFLookup is empty — provide a density reduction factor table in params.")
  # terra::values<- deep-copies before writing, so studyAreaRaster is not modified
  sim$snagBiomass_Mg_ha <- sim$studyAreaRaster
  terra::values(sim$snagBiomass_Mg_ha) <- NA_real_
  sim$DWDBiomass_Mg_ha  <- sim$studyAreaRaster
  terra::values(sim$DWDBiomass_Mg_ha)  <- NA_real_
  sim$snagHistory <- NULL
  sim$DWDHistory  <- NULL
  return(invisible(sim))
}

Annual <- function(sim) {
  drf <- P(sim)$DRFLookup

  # Right join: inventory rows with no matching DRF get DRF=NA; na.rm=TRUE silently drops them
  snagWithBiomass <- drf[pool == "snag"][sim$snagTable, on = .(species, DC)]
  snagWithBiomass[, currentBiomass := initBiomass * DRF]
  snagByPixel <- snagWithBiomass[, .(value = sum(currentBiomass, na.rm = TRUE)), by = pixelID]
  sim$snagBiomass_Mg_ha <- pixelValuesToRaster(snagByPixel, sim$studyAreaRaster)

  DWDwithBiomass <- drf[pool == "DWD"][sim$DWDTable, on = .(species, DC)]
  DWDwithBiomass[, currentBiomass := initBiomass * DRF]
  DWDbyPixel <- DWDwithBiomass[, .(value = sum(currentBiomass, na.rm = TRUE)), by = pixelID]
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
  return(invisible(sim))
}
