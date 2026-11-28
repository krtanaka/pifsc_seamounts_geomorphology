# PIFSC Seamounts & Geomorphology

Distance-to-seamount rasters and seafloor geomorphology layers for use as
predictors in Pacific longline fisheries models (CPUE, catch composition,
bycatch).

| Script | Purpose |
|---|---|
| `process_knolls_seamounts.R` | Filters Yesson et al. (2011) seamounts/knolls and builds global geodesic distance rasters (NetCDF) + maps |
| `process_geomorphology_shapefiles.R` | Rasterizes Harris et al. (2014) global seafloor geomorphology |

## Simple raster-only scripts

These standalone scripts need only `terra` (`install.packages("terra")`).
They use the staged dataset folders under `~/data/deepsea/data_download`.
Set the environment variable `PIFSC_GIS_DATA_DIR` to use another root directory,
or edit `data_root` at the top of each script.

| Script | Output |
|---|---|
| `rasterize_seamounts_knolls.R` | Distance in km to the nearest summit; edit `max_depth_m` to `400`, `1000`, or `NULL` (all depths) |
| `rasterize_geomorphology.R` | One estimated fractional-coverage GeoTIFF (0-1) per geomorphology shapefile |
| `rasterize_ecological_benthic_units.R` | Existing EBU GeoTIFFs aligned to the common grid, without re-rasterizing duplicate vector products |

All use WGS84, longitude -180 to 180, latitude -60 to 60, at 0.05 degrees
(7,200 columns by 2,400 rows). Outputs go under `output/rasters/`, separated by
dataset. Existing outputs are skipped; set `overwrite <- TRUE` after changing
inputs or options. The scripts run sequentially, without maps or worker clusters.
Seamount distances use all global input points to avoid boundary effects; the
400 m and 1,000 m thresholds are inclusive. Surface-breaking summits are removed;
missing depths are retained only for the all-depth run. Confirm the inferred
`DEPTH` sign and metre units against source metadata.

From the repository directory on the cloud workstation:

```bash
Rscript rasterize_seamounts_knolls.R
Rscript rasterize_geomorphology.R
Rscript rasterize_ecological_benthic_units.R
```

EBU reads only `01_classification`, `02_environmental_predictors`, and
`04_supporting`; it excludes BioOracle, temporary/smoke-test folders and partial
downloads. Classification, lithology and slope-hemisphere rasters use
nearest-neighbor resampling; continuous rasters use bilinear interpolation.
Source units are unchanged. Reprojection does not make every source suitable
for modeling: the EBU metadata reports a suspicious 0-360 range for
`EBU_predictor_bathymetry_0bd2bd53.tif`; the script warns but preserves that layer
under its original name for review. Several predictors have multiple versions;
these remain separate, not silently merged.

---

## 1. Should seamounts be filtered? Literature notes

### Key evidence

- **Morato et al. 2010, *PLOS ONE*** — "Tuna longline fishing around West and
  Central Pacific seamounts." doi:[10.1371/journal.pone.0014453](https://doi.org/10.1371/journal.pone.0014453)
  - Seamount definition: features >= 100 m high that do not break the surface
    (Kitchingman & Lai database, 7,741 features, validated by Allain et al.).
  - Analysed longline sets within 100 km of a seamount; estimated seamount
    influence at **~20–40 km**.
  - Effects are not uniform: 36% of screened seamounts showed higher catch rates
    for at least one tuna species, 41% showed *lower* rates.
  - Summit depth could not be modelled (data limitations), but the authors note
    many deep seamounts are unlikely to aggregate pelagic species.
  - A single longline set can extend > 100 km, which limits spatial precision.
- **Morato et al. 2010, *PNAS*** — doi:[10.1073/pnas.0910290107](https://doi.org/10.1073/pnas.0910290107)
  - Pelagic diversity higher within **30–40 km** of summits.
  - Higher catches near seamounts: yellowfin tuna, swordfish, blue marlin,
    silky, mako and porbeagle sharks. Lower: albacore, shortbill spearfish.
- **Lima et al. 2020, *Front. Mar. Sci.*** — doi:[10.3389/fmars.2020.00238](https://doi.org/10.3389/fmars.2020.00238)
  - NE Atlantic mechanistic context; small seamounts 100–1,000 m high, large
    > 1,000 m. Productivity/aggregation effects linked to summit shallowness.
- **Gilman et al. 2012, *Aquatic Conserv.*** — doi:[10.1002/aqc.2237](https://doi.org/10.1002/aqc.2237)
  (seamount–longline interactions; not fully reviewed, access blocked).
- **Wessel et al. 2010, *Oceanography*** — doi:[10.5670/oceanog.2010.60](https://doi.org/10.5670/oceanog.2010.60)
  (seamount census, definitions).
- **Yesson et al. 2011, *Deep-Sea Res. I*** — doi:[10.1016/j.dsr.2011.02.004](https://doi.org/10.1016/j.dsr.2011.02.004)
  (source dataset used here; ~33k seamounts > 1,000 m high and ~138k knolls
  200–1,000 m high).

### Is "summit < 400 m" reasonable?

**Yes, as the primary biologically motivated layer — but not as a universal
threshold.** Rationale:

- Mechanisms that aggregate pelagics (Taylor columns, upwelling, enhanced
  primary production, trapping of vertically migrating prey) weaken with summit
  depth; shallow summits interact with the epipelagic/upper mesopelagic layer
  where longline hooks fish (shallow-set swordfish ~30–90 m; deep-set tuna
  ~100–400 m).
- No study we found tests 400 m as a cutoff explicitly; Morato et al. could not
  model summit depth. So treat 400 m as a hypothesis and test sensitivity.

### Recommended layers

| Layer | Filter | Role |
|---|---|---|
| **Primary** | summit depth <= 400 m, seamounts + knolls | main predictor |
| **Sensitivity** | summit depth <= 1,000 m | robustness check (our proposal, not a published cutoff) |
| **Baseline** | all features | comparison / null |

### Other recommendations

- Filter on **summit depth (`DEPTH`)**, not `HEIGHT`. Height is relief above the
  surrounding seafloor; a 3,000 m tall seamount can still have a 2,500 m deep
  summit. The previous `HEIGHT >= 2500` filter is far stricter than the 100 m
  definition and can exclude shallow features.
- **Include knolls** initially: Yesson et al. knolls (200–1,000 m relief) can
  have shallow summits relevant to pelagics.
- Check the sign/units of `DEPTH` (the script infers it; see options).
- Drop features with missing depth (when depth-filtering), summits above sea
  level (islands/atolls), and duplicate locations.
- Use **geodesic distance to the summit point**; include features outside the
  study area to avoid edge effects (the script computes globally).
- The 5 km grid is `5/111.12` degrees (~0.045°), so cells are ~5 km only at the
  equator; east–west width shrinks with latitude.
- Model by species and fishery (deep-set vs shallow-set), allow non-linear
  effects, and use spatially blocked cross-validation.
- Treat 30–40 km as an exploratory scale (e.g. for a binary "near seamount"
  variable), not a hard cutoff; distance itself is the preferred predictor.

---

## 2. Running `process_knolls_seamounts.R`

Edit only the **OPTIONS** block at the top, then source the script.

### Inputs

`base_path` must point to the Yesson et al. (2011) `ModelledSeamounts2011_v1/01_Data/`
folder containing `Seamounts/Seamounts.shp` and `Knolls/Knolls.shp`
(fields used: `HEIGHT`, `DEPTH`).

### Filter scenarios

Each entry in `scenarios` produces its own NetCDF(s). Unset fields use
`scenario_defaults`.

| Option | Meaning | Default |
|---|---|---|
| `include_seamounts` | use the Seamounts layer | `TRUE` |
| `include_knolls` | use the Knolls layer | `TRUE` |
| `max_summit_depth_m` | keep summit depth <= value (m, positive down); `NULL` = no filter | `NULL` |
| `min_height_m` | keep `HEIGHT` >= value (m); `NULL` = no filter | `NULL` |
| `height_tertiles` | also write Small/Medium/Big rasters (HEIGHT tertiles within the filtered set) | `FALSE` |

Default scenarios: `primary_shallow` (<= 400 m), `sensitivity_1000m`
(<= 1,000 m), `baseline_all`. A commented-out `legacy_tall` scenario
reproduces the previous version (seamounts only, `HEIGHT >= 2500`, tertiles).

### Other options

| Option | Meaning |
|---|---|
| `depth_sign_convention` | `"auto"` (infer from median), `"negative_down"`, `"positive_down"` |
| `drop_duplicate_locations` | drop identical summit coordinates |
| `recalculate` | `"missing"` compute only absent NetCDFs (default), `"all"` overwrite, `"none"` plot only |
| `res_km`, `n_chunks`, `n_cores` | grid resolution and parallelisation; workers automatically respect CPU quotas, reserving two cores (minimum one worker) and capped by the number of strips |
| `make_feature_maps` | exploratory maps/hex plots of raw layers (hex plots mark depth thresholds) |
| `make_distance_maps` | one log10 distance map per NetCDF |
| `nc_dir`, `output_dir` | output locations (defaults: repo root, `output/`) |

### Output naming

```
Distance_<Features>_<Depth>_<Height>_<SizeClass>_<Res>km.nc
Map_<Features>_<Depth>_<Height>_<SizeClass>_<Res>km_Log10.png
```

- `Features`: `Seamounts`, `Knolls`, or `SeamountsKnolls`
- `Depth`: `DepthLE400m` (summit <= 400 m) or `DepthAll`
- `Height`: `HgtGE2500m` (height >= 2500 m) or `HgtAll`
- `SizeClass`: `All`, `Small`, `Medium`, `Big`

Examples:

- `Distance_SeamountsKnolls_DepthLE400m_HgtAll_All_5km.nc` (primary)
- `Distance_SeamountsKnolls_DepthLE1000m_HgtAll_All_5km.nc` (sensitivity)
- `Distance_SeamountsKnolls_DepthAll_HgtAll_All_5km.nc` (baseline)
- `Distance_Seamounts_DepthAll_HgtGE2500m_Big_5km.nc` (legacy)

Each NetCDF stores variable `distance` (km, 0–360 longitudes) with the filter
settings, feature count and source citation as global attributes.
`output/distance_raster_manifest.csv` lists every scenario, its filters,
feature counts and file path.

The older files `Distance_{All,Small,Medium,Big}_Seamounts_5km.nc` were produced
by the previous version (seamounts with `HEIGHT >= 2500`, height tertiles) and
correspond to the `legacy_tall` scenario.

### Runtime note

Both scripts require `parallelly` for cloud/container-aware CPU detection.
Install it once in R with `install.packages("parallelly")`. The geomorphology
script likewise reserves two cores and caps workers at the number of layers.
Automatic detection replaces the previous 20-worker cap and fixed five-worker
override, respectively. CPU detection does not account for raster memory needs;
lower `n_cores` or `num_cores` manually if memory is limited.

Global 5 km distance rasters are expensive (minutes to hours per layer,
depending on cores and feature count). Leave `recalculate = "missing"` so
re-runs only compute new scenarios.
