# Future implementations

Ideas that have been scoped and deliberately deferred. Each entry says what it is, why it
was deferred, and what it depends on, so it can be picked up without re-deriving the context.

## Environmental forcing along real routes (biofouling model)

**Status:** deferred on 2026-11-07, when the global-movement redesign of the fouling model
was scoped to the model only ("Option A"). This entry is "Option B".

**What:** build the pipeline that turns a real ship track into the `Segment`s the fouling
model consumes, so every segment carries the environment the ship actually experienced:

- **Track → segments.** Split AIS positions into segments of near-constant conditions
  (speed through water, port status, distance to coast), with a maximum segment length of
  about one day so the model's segment-end observations stay accurate (see "Segmentation"
  in `docs/methodology/ocean-fouling.md`). Detect berth and anchorage stays (Global Fishing
  Watch anchorages) and dry-dock visits.
- **Sample ocean fields** at each segment's time and position:

  | Variable | Source |
  |---|---|
  | Sea-surface temperature | NOAA OISST v2.1 (daily, 0.25°) |
  | Sea-surface salinity | Copernicus GLORYS12 |
  | Chlorophyll-a | Copernicus GlobColour L4 (daily or monthly, 4 km) |
  | Distance to coast | NASA/GSFC distance-to-coast grid |

  Use time-and-space interpolation with a defined fallback where a value is missing (sea ice,
  cloud gaps in chlorophyll).
- **Speed through water** from speed over ground and ocean currents (Copernicus currents),
  because settlement and friction depend on flow past the hull, not on speed over ground.
- **Storage and caching** of the sampled fields, so re-running a fleet does not re-download
  them. Everything must be deterministic and reproducible from pinned dataset versions.

**Why deferred:** it is a data-engineering project of its own (downloads, storage,
interpolation, AIS cleaning). The model has to handle arbitrary global forcing correctly
first, and synthetic world routes are enough to test that.

**Depends on:** the global-movement fouling model (Option A). It also needs the AIS ingestion
in the ocean pipeline (plan 0001) and a decision on local vs remote storage for gridded
datasets. All of the sources above are free, in line with the project's no-licensing-fees rule.

**Spec references:** `research/reports/Biofouling growth dynamical model.md`, sections on
environmental drivers and data sources (the "Anchorages, port visits" and "Chlorophyll-a"
rows of the data table).
