DeadWood_Biomass
================
2026-05-11

- [DeadWood_Biomass](#deadwood_biomass)
  - [Inputs](#inputs)
  - [Outputs](#outputs)
  - [Parameters](#parameters)
  - [Events](#events)
    - [`init` (once, at simulation
      start)](#init-once-at-simulation-start)
    - [`transition` (every 5 years, priority
      4)](#transition-every-5-years-priority-4)
    - [`plot` (optional, every `.plotInterval` years, priority
      5)](#plot-optional-every-plotinterval-years-priority-5)
  - [Event scheduling and module
    interactions](#event-scheduling-and-module-interactions)
  - [Usage example](#usage-example)
  - [Package dependencies](#package-dependencies)

# DeadWood_Biomass

A [SpaDES](https://spades.predictiveecology.org/) module that converts
the snag and downed woody debris (DWD) decay class inventories into
pixel-level biomass estimates (Mg ha⁻¹). It reads the current
`snagTable` and `DWDTable` produced by `DeadWood_snagDecay` and
`DeadWood_DWDDecay`, aggregates biomass by pixel, and writes the results
to spatial rasters. Optionally, it accumulates a time-series of biomass
snapshots for visualization and analysis.

------------------------------------------------------------------------

## Inputs

| Object | Class | Description |
|----|----|----|
| `snagTable` | `data.table` | Current snag inventory from `DeadWood_snagDecay`. Must have columns `pixelID` (integer) and `initBiomass` (numeric, Mg ha⁻¹). |
| `DWDTable` | `data.table` | Current DWD inventory from `DeadWood_DWDDecay`. Must have columns `pixelID` (integer) and `initBiomass` (numeric, Mg ha⁻¹). |
| `studyAreaRaster` | `SpatRaster` | Template raster defining the spatial grid (extent, CRS, resolution). All output rasters are created to match this template. |

------------------------------------------------------------------------

## Outputs

| Object | Class | Description |
|----|----|----|
| `snagBiomass_Mg_ha` | `SpatRaster` | Pixel-level snag current biomass (Mg ha⁻¹), updated every 5 years. Integrates DRF by decay class. Pixels with no snags have value `NA`. |
| `DWDBiomass_Mg_ha` | `SpatRaster` | Pixel-level DWD current biomass (Mg ha⁻¹), updated every 5 years. Integrates DRF by decay class. Pixels with no DWD have value `NA`. |
| `snagHistory` | `SpatRaster` | Multi-layer raster accumulating one snag biomass snapshot per plot event (layers named `yr<time>`). `NULL` if `.plotInitialTime` is `NA`. |
| `DWDHistory` | `SpatRaster` | Multi-layer raster accumulating one DWD biomass snapshot per plot event (layers named `yr<time>`). `NULL` if `.plotInitialTime` is `NA`. |

------------------------------------------------------------------------

## Parameters

| Parameter | Type | Default | Description |
|----|----|----|----|
| `.plotInitialTime` | numeric | `NA` | Simulation time of the first snapshot. Set to `NA` to disable snapshot accumulation entirely. |
| `.plotInterval` | numeric | `5` | Interval (years) between snapshots. Defaults to 5 to align with the decay timestep. Only used if `.plotInitialTime` is not `NA`. |
| `DRFLookup` | `data.table` | *Pinus strobus* defaults | Density reduction factors by species, pool (`"snag"` or `"DWD"`), and DC. Columns: `species` (chr), `pool` (chr), `DC` (int), `DRF` (num). Default values from Paper 2 Appendix D. |

------------------------------------------------------------------------

## Events

### `init` (once, at simulation start)

- Creates `snagBiomass_Mg_ha` and `DWDBiomass_Mg_ha` as copies of
  `studyAreaRaster` filled with `NA`.
- Initialises `snagHistory` and `DWDHistory` to `NULL`.
- Schedules the first `transition` event at `start(sim) + 5` with
  priority 4.
- If `.plotInitialTime` is not `NA`, schedules the first `plot` event at
  `.plotInitialTime` with priority 5.

### `transition` (every 5 years, priority 4)

Runs after `DeadWood_snagDecay` and `DeadWood_DWDDecay` have updated
their inventories within the same timestep. Computes pixel-level biomass
sums for both pools.

**Biomass calculation:**

For each pool, `initBiomass` (biomass at time of death) is multiplied by
the decay-class-specific density reduction factor (DRF) and summed by
pixel:

$$\text{snagBiomass}[p] = \sum_{i \,:\, \text{pixelID}_i = p} \text{initBiomass}_i \times \text{DRF}[\text{species}_i,\; \text{snag},\; \text{DC}_i]$$

$$\text{DWDBiomass}[p] = \sum_{i \,:\, \text{pixelID}_i = p} \text{initBiomass}_i \times \text{DRF}[\text{species}_i,\; \text{DWD},\; \text{DC}_i]$$

where $p$ is a pixel index, $i$ indexes individual pieces, and DRF
values come from the `DRFLookup` parameter. DRF \< 1 for DC \> 1
reflects the reduction in wood density as decomposition progresses.
Pixels with no pieces receive `NA`.

The resulting pixel vectors are mapped onto a copy of `studyAreaRaster`
using `pixelValuesToRaster()`, which sets all pixels not present in the
table to `NA`.

### `plot` (optional, every `.plotInterval` years, priority 5)

Captures snapshots of `snagBiomass_Mg_ha` and `DWDBiomass_Mg_ha` at the
current simulation time and appends them as named layers to
`snagHistory` and `DWDHistory`. Layer names follow the pattern
`yr<time>` (e.g., `yr10`, `yr25`).

Snapshot events continue until the end of the simulation. The last
snapshot scheduled after `end(sim)` is not fired.

------------------------------------------------------------------------

## Event scheduling and module interactions

This module runs at the lowest priority within each 5-year timestep,
after both decay modules have updated their inventories:

| Priority | Module | Event | Purpose |
|----|----|----|----|
| 1 | `DeadWood_snagDecay` | `transition` | Advance snag DC, produce `fallenSnags` |
| 2 | `DeadWood_DWDDecay` | `receive` | Accept `fallenSnags` into DWD pool |
| 3 | `DeadWood_DWDDecay` | `transition` | Advance DWD DC via logistic model |
| 4 | `DeadWood_Biomass` | `transition` | Compute biomass rasters ← **this module** |
| 5 | `DeadWood_Biomass` | `plot` | Capture biomass snapshots ← **this module** |

------------------------------------------------------------------------

## Usage example

``` r
library(SpaDES.core)
library(terra)

# Minimal 3x3 study area raster
studyAreaRaster <- terra::rast(
  nrows = 3, ncols = 3,
  xmin = 0, xmax = 3, ymin = 0, ymax = 3
)

# snagTable and DWDTable would normally come from the decay modules
snagTable <- data.table::data.table(
  pixelID     = c(1L, 1L, 5L),
  species     = "Pinus strobus",
  DC          = c(1L, 2L, 3L),
  ageInDC     = c(0L, 5L, 10L),
  initBiomass = c(12.5, 8.3, 6.0),
  diameter_cm = c(22.0, 18.5, 15.0)
)

DWDTable <- data.table::data.table(
  pixelID       = c(2L, 5L),
  species       = "Pinus strobus",
  DC            = c(3L, 2L),
  ageInDC       = c(0L, 5L),
  ageSinceEntry = c(0L, 5L),
  initBiomass   = c(10.0, 7.5),
  diameter_cm   = c(20.0, 16.0)
)

mySim <- simInit(
  times   = list(start = 0, end = 50),
  params  = list(DeadWood_Biomass = list(
    .plotInitialTime = 10,
    .plotInterval    = 10
  )),
  modules = list("DeadWood_Biomass"),
  objects = list(
    snagTable        = snagTable,
    DWDTable         = DWDTable,
    studyAreaRaster  = studyAreaRaster
  )
)

mySim <- spades(mySim)

# Current biomass rasters
mySim$snagBiomass_Mg_ha
mySim$DWDBiomass_Mg_ha

# Time-series snapshots (if .plotInitialTime was set)
mySim$snagHistory
terra::plot(mySim$snagHistory)
```

------------------------------------------------------------------------

## Package dependencies

- [`SpaDES.core`](https://github.com/PredictiveEcology/SpaDES.core) (\>=
  3.0.0)
- [`data.table`](https://CRAN.R-project.org/package=data.table)
- [`terra`](https://CRAN.R-project.org/package=terra)
- [`ggplot2`](https://CRAN.R-project.org/package=ggplot2)
