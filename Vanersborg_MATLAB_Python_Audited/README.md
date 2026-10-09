# Vänersborg MATLAB and Python pipelines — audited recovery

Keep all three MATLAB files together:
- `vanersborg_full_pipeline.m`: compatibility entry point.
- `vanersborg_full_pipeline_from_zip.m`: complete MATLAB workflow, tables and figures.
- `vanersborg_core.m`: shared numerical routines executed in GNU Octave 8.4.0.

The Python workflow is `vanersborg_full_pipeline_from_zip.py`.

## MATLAB

Extract into a new folder and make that folder MATLAB's current folder.
Use MATLAB R2021b or later. No Signal Processing Toolbox is required.

```matlab
clear vanersborg_full_pipeline vanersborg_full_pipeline_from_zip vanersborg_core
rehash
which vanersborg_full_pipeline -all
which vanersborg_full_pipeline_from_zip -all
which vanersborg_core -all
vanersborg_full_pipeline( ...
    'D:\data\Dataset_steel_bridge_in_Sweden_1.zip', ...
    'D:\data\Dataset_steel_bridge_in_Sweden_2.zip', ...
    'D:\results\Vanersborg_audited');
```

Alternatively, put both raw ZIP archives alongside the MATLAB files and run
`vanersborg_full_pipeline` without arguments. Use a new output folder to avoid
mixing earlier and corrected results. The raw archives are not bundled here.

## Python

```bash
python -m pip install -r requirements.txt
python vanersborg_full_pipeline_from_zip.py --zip1 /data/Dataset_steel_bridge_in_Sweden_1.zip --zip2 /data/Dataset_steel_bridge_in_Sweden_2.zip --out results --workers 4
```

After a completed alignment, `--reuse-cache` reuses the aligned NPZ in the output
folder. `--reuse-extracted` reuses extracted CSVs.

## What was executed

The Python workflow was run from both original raw archives (64 events), including
tables and PNG/PDF figures. The MATLAB numerical core was executed in Octave on
all 64 raw numeric records. It was compared against Python for alignment, features,
POD/QR selection, virtual strains, event scores, primary counts and all 20 baseline
partitions. See `AUDIT_REPORT.md` and `audit/` for measured differences and logs.

The MATLAB `readtable`/`table`/`writetable` and graphics/export layer was not run:
MATLAB is unavailable and Octave does not implement all these MATLAB APIs.
This is an executed numerical-core audit, not a claim of a full native MATLAB run.

## Reproduce the Octave numerical audit

First run the Python pipeline above. Fixture preparation reads its extracted raw
CSVs and stores the selected raw numeric columns without alignment; Octave itself
performs alignment, fitting and inference.

```bash
python prepare_octave_audit.py --python-results results --out audit_fixtures
octave --no-gui --quiet --eval "addpath(pwd); audit_vanersborg_octave('audit_fixtures','octave_results');"
python compare_octave_audit.py --python-results results --octave-results octave_results
```

The full MATLAB workflow and audit both call the same `vanersborg_core.m` functions.
The audit uses Octave's MAT-file save option and is intended to run in Octave.
Bootstrap draws use each language's own random generator, so individual resampled
bootstrap outcomes need not be identical across languages. Primary record counts
and all 20 exhaustive partitions are compared exactly.

## Results and interpretation

| Configuration | Detected | Correct top zone | Detected and correct |
|---|---:|---:|---:|
| 12 physical gauges | 7/7 | 7/7 | 7/7 |
| 5 QR physical gauges | 2/7 | 7/7 | 2/7 |
| 5 QR + 7 virtual gauges | 7/7 | 7/7 | 7/7 |

QR references: SG11, SG3, SG1, SG13, SG10. SG9 is reconstructed; SG10 is measured.
Primary fitting/calibration uses the first three/next three early events. Virtual
strains are healthy-trained reconstructions, not independent measurements. The
seven crack/post-crack records and 20 overlapping baseline partitions do not
constitute independent validation trials.

`verified_python_run/` contains newly generated outputs, not inherited old results.
