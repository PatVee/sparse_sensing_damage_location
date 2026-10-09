# Sparse strain sensing for bridge damage localisation

Reproducible scripts for two complementary studies: **KW51** (numerical stiffness-loss localisation using measured healthy train responses and an FE model) and **Vänersborg** (experimental crack-region assessment using physical, QR-reduced and reconstructed virtual strain sensors).

## KW51 — MATLAB

Work in the folder containing `RUN_KW51.m` (`KW51_MATLAB_PACKAGE/KW51_MATLAB_PACKAGE/` after extraction).

1. Use **MATLAB R2022b+**, with Signal Processing Toolbox and the MATLAB JVM.
2. Place `traindata_201901_part1.zip` and `traindata_201901_part2.zip` in `input/`, or select them when prompted. Leave the ZIPs compressed; the FE reference data are already in `reference/`.
3. Run in MATLAB:

   ```matlab
   RUN_KW51
   ```

The script extracts and aligns healthy passages, builds FE damage signatures, evaluates sparse strain features and reports detection/localisation across candidate zones. Adjust settings in `KW51_CONFIG.m`. Results appear in `outputs/tables/`, `outputs/trials/`, `outputs/cache/` and `outputs/figures/`.

To **redraw figures without rerunning the analysis**, use `PRINT_KW51_RESULTS`. The separate `full_make_figure.m` also works from saved output tables/trials.

## Vänersborg — MATLAB or Python

Work in the extracted `Vanersborg_MATLAB_Python_Audited/` folder. Obtain the two original raw archives: `Dataset_steel_bridge_in_Sweden_1.zip` and `Dataset_steel_bridge_in_Sweden_2.zip` (not included).

**MATLAB (R2021b+):** keep `vanersborg_full_pipeline.m`, `vanersborg_full_pipeline_from_zip.m` and `vanersborg_core.m` together, then run:

```matlab
vanersborg_full_pipeline('Dataset_steel_bridge_in_Sweden_1.zip', ...
    'Dataset_steel_bridge_in_Sweden_2.zip', 'results')
```

**Python:**

```bash
python -m pip install -r requirements.txt
python vanersborg_full_pipeline_from_zip.py --zip1 Dataset_steel_bridge_in_Sweden_1.zip --zip2 Dataset_steel_bridge_in_Sweden_2.zip --out results
```

The workflow compares the full physical array, **five QR-selected physical gauges**, and those same references augmented by **seven healthy-trained virtual strain reconstructions**. It produces CSV summaries and PDF/PNG plots in the selected results folder. Virtual gauges are reconstructed estimates, **not independent measurements or FE-generated sensors**.

## Interpretation and reproducibility

Both studies evaluate spatial strain information using area/Gram/DCT-based features, but their evidence is different: **KW51 damage cases are FE-prescribed**, while **Vänersborg uses an inspection-confirmed crack region**. The KW51 package documents code/dependency checks but not a completed native MATLAB numerical rerun; the Vänersborg package includes an executed Python run and an Octave numerical-core comparison (see its `AUDIT_REPORT.md`), not a full native MATLAB execution.

**Note:** KW51 `RUN_KW51.m` is the *sparse-localisation* entry point; it does **not** require the separate probabilistic virtual-sensing `main.m` or `main_kw51_modal_only.m` workflow.
