# KW51 MATLAB analysis package

## Start

1. Extract this entire ZIP, keeping its folders together.
2. Put `traindata_201901_part1.zip` and `traindata_201901_part2.zip` in `input/`,
   or select their existing locations when prompted.
3. Open `RUN_KW51.m` in MATLAB and press **Run**.

Only the two train archives are supplied by you. The FE model, candidate
geometry and sensor coordinates are already included in `reference/`; the
main runner locates them relative to itself. No `D:\model` paths or Python
installation are needed.

MATLAB R2022b or later, Signal Processing Toolbox and the MATLAB JVM are required.
Change trial counts, passage windows and output location in `KW51_CONFIG.m`.

## File responsibilities

| File | Responsibility |
|---|---|
| `RUN_KW51.m` | Main runner: select archives, call analysis, then print results |
| `KW51_CONFIG.m` | Input/output locations and analysis settings |
| `src/kw51_run_analysis.m` | Coordinate stages and validated checkpoints |
| `src/kw51_extract_traindata.m` | Extract ZIPs, deduplicate, screen channels, resample, remove baseline, filter and align passages |
| `src/kw51_build_fe_library.m` | Full-matrix static strain calculation for healthy and damaged physical zones |
| `src/kw51_defineZones.m` | Construct the frozen 68-zone shortlist from bundled FE geometry |
| `src/kw51_build_features.m` | Build measured features, damage signatures and QR selections |
| `src/kw51_featureOne.m` | Explicit feature definitions: RMS/P2P, signed area, Gram, DCT and combinations |
| `src/kw51_fitHealthy.m` | Fit healthy class/temperature mean and regularised covariance |
| `src/kw51_localisation_trials.m` | Calibrate detection threshold and simulate repeated-passage localisation trials |
| `src/kw51_write_tables.m` | Trial counts, probabilities, Wilson intervals, feature rankings and retained locations |
| `PRINT_KW51_RESULTS.m` | Read generated tables/trial files and export all plots |
| Other `src/*.m` files | FE beam, loading, feature and file utilities called by these stages |

## Why the model files are included

The train archives supply **measured healthy response variability**. They do
not supply the bridge stiffness matrix, its mesh connectivity or the strain
operator needed to simulate damage. `reference/base_and_all_gauge_calibration.mat`
contains the previously supplied assembled KW51 export; the solver uses its
`Data.Base` fields. `candidate_regions.csv` contains geometric/model metadata,
not previously computed localisation results. `strain_sensor_locations.csv`
contains the positions for the 11 active measured strain channels. Paired
gauges sharing an FE operator-node centroid overlap in the red-dot plot.

No precomputed damage library or localisation-result CSV is required as an input.
The bundled reference CSVs are model metadata and do not replace raw extraction.

## Output folders

- `outputs/tables/`: inventories, allocations, channel quality, geometry,
  feature catalogue/QR indices, observability diagnostics, GLRT results,
  probabilities, counts, Wilson intervals and rankings.
- `outputs/tables/per_location/`: every feature and severity for each tested zone.
- `outputs/trials/`: individual predictions, detections, scores and resampled
  event IDs for every feature/window. Calibration draws are also saved.
- `outputs/cache/`: aligned responses, FE passages, feature libraries and GLRT
  results. Input/source/settings fingerprints prevent reuse of mismatched caches.
- `outputs/figures/`: PDF and PNG retained-zone overview, localisation
  probability comparison and trial-outcome bars.
- `outputs/figures/per_location/`: detection, ungated localisation and
  detection-gated localisation curves with Wilson intervals for all 68 zones
  by default. The feature is selected for each complete 10--50% curve.
- `outputs/figures/spatial_predictions/`: retained-zone 3D predictions for
  every severity and passage window, with sensor locations shown as red dots.
- `outputs/work/`: extracted archives and intermediate processed records.

Run `PRINT_KW51_RESULTS` on its own to redraw existing output. It reads only
generated tables/trial matrices; it does not extract traindata or run the solver.
For an alternative output folder, use `PRINT_KW51_RESULTS('D:\my_results')`.

All 13 feature families and 68 zones are kept in the tables. Retention is based
on a minimum detected-and-correct fraction of 0.50 across all five severities
in either default window, with **no cap on the number of retained locations**.
All tied best features by state and all Pareto features are separately exported.

## Protocol and interpretation

The fixed chronological split is 32 healthy fitting records, 14 threshold
calibration records and 16 later residual records. The function checks for the
original 62 unique January events and refuses to silently apply that split to
a different dataset. Eleven strain channels are used; `sgBD1718C` is excluded
because it was unavailable in this dataset. This final location-search protocol
does not use acceleration features. Filtering is 2 Hz after 50 Hz resampling.

The nominal 12-axle train is applied on both tracks in the static FE model.
Local stiffness loss is 10%, 20%, 30%, 40% or 50%, applied to physical zones.
This rebuilds the frozen 68-zone shortlist used in the final feature search;
it does not repeat the earlier exhaustive 1,603-element screening study.
The response-class rule is an RMS-ratio surrogate, not measured track metadata.

The default 300 trials per state bootstrap 16 or 32 passages from the 16
later residual records. They are not 300 independent measured damaged events.
Truth and inverse signatures use the same FE model. The probabilities are
empirical performance fractions for this exploratory semi-synthetic experiment.
Wilson intervals describe conditional Monte Carlo sampling variability, not
field reliability or uncertainty in the physical damage model. The 95th
percentile calibration threshold does not establish 95% operational specificity.
Location and feature selection on these trials is exploratory, not blind validation.

In the 3D prediction figures, bold segments identify the prescribed zone and
open-circle area represents the fraction of **all trials** that detected and
predicted each zone. Predictions are located at FE zone centroids. Red dots
are available strain sensor positions, not inferred damage coordinates.

## Validation status

The MATLAB source was syntax/lint checked and the package dependencies were
checked. No traindata archives were opened or processed while preparing this
package, and the MATLAB numerical workflow has not been executed here.
The previously reported CSVs were generated in Python. MATLAB's random stream,
filtering numerics and QR tie ordering can cause differences in a fresh run;
this package does not claim bit-for-bit reproduction of earlier trial counts.
