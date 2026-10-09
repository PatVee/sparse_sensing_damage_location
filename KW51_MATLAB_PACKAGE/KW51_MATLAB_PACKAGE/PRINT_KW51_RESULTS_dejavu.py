"""
Python equivalent of PRINT_KW51_RESULTS.m using DejaVu Sans.

The script redraws KW51 result figures from saved tables and saved trial
matrices. It does NOT reopen the raw train data and does NOT rerun the FE
solver or Monte Carlo trials.

Expected structure
------------------
<output_root>/
    tables/
        candidate_regions.csv
        zone_definitions.csv
        retained_best_locations.csv
        best_feature_trials_for_every_location.csv
        strain_sensor_locations.csv
    trials/
        <feature>_<n>pass.mat
    figures/                  # created automatically
        per_location/
        spatial_predictions/

Examples
--------
python PRINT_KW51_RESULTS_dejavu.py
python PRINT_KW51_RESULTS_dejavu.py --output-root "D:/results/KW51/outputs"
python PRINT_KW51_RESULTS_dejavu.py --output-root "D:/results/KW51/outputs" --retained-only
"""

from __future__ import annotations

import argparse
import io
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# -----------------------------------------------------------------------------
# Font and publication style
# Same font family used in the Python truss plotting script.
# -----------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.18,
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})


def gray(value: float):
    return (value, value, value)


def resolve_output_root(user_root: str | None) -> Path:
    if user_root:
        root = Path(user_root).expanduser().resolve()
    else:
        # MATLAB default:
        # outputRoot = fullfile(fileparts(mfilename('fullpath')),'outputs')
        root = Path(__file__).resolve().parent / "outputs"

    tables = root / "tables"
    if not tables.is_dir():
        raise FileNotFoundError(
            "Could not find the expected tables directory:\n"
            f"  {tables}\n"
            "Use --output-root to specify the executed KW51 output folder."
        )
    return root


def save_figure(fig, path_without_suffix: Path, dpi: int = 220):
    """
    Robust PDF/PNG export using an in-memory buffer before writing to disk.
    """
    path_without_suffix.parent.mkdir(parents=True, exist_ok=True)

    for ext in ("png", "pdf"):
        stream = io.BytesIO()
        fig.savefig(
            stream,
            format=ext,
            dpi=dpi,
            bbox_inches="tight",
            pad_inches=0.04,
            facecolor="white",
        )
        payload = stream.getvalue()

        if ext == "pdf":
            if not (payload.startswith(b"%PDF") and b"%%EOF" in payload[-200:]):
                raise RuntimeError(f"PDF export validation failed for {path_without_suffix}.pdf")

        (path_without_suffix.parent / f"{path_without_suffix.name}.{ext}").write_bytes(payload)


def read_tables(output_root: Path):
    tables = output_root / "tables"

    regions = pd.read_csv(tables / "candidate_regions.csv")
    zones = pd.read_csv(tables / "zone_definitions.csv")
    retained = pd.read_csv(tables / "retained_best_locations.csv")
    trials = pd.read_csv(tables / "best_feature_trials_for_every_location.csv")
    sensors = pd.read_csv(tables / "strain_sensor_locations.csv")

    required_region = {"elementID", "x1_m", "y1_m", "z1_m", "x2_m", "y2_m", "z2_m"}
    required_zone = {"zone_index", "zone_name", "element_ids", "x_m", "y_m", "z_m"}
    required_retained = {"zone_index", "zone_name", "element_ids", "x_m", "y_m", "z_m"}
    required_sensor = {"x_m", "y_m", "z_m"}
    required_trials = {
        "zone_index", "n_passages", "severity_pct", "feature", "loss",
        "detection_rate", "localisation_rate",
        "detected_and_correctly_localised_rate",
        "detected_trials", "correctly_localised_trials",
        "detected_and_correctly_localised_trials",
        "detected_wrong_location_trials",
        "correctly_localised_without_detection_trials",
        "total_trials",
    }

    checks = [
        ("candidate_regions.csv", regions, required_region),
        ("zone_definitions.csv", zones, required_zone),
        ("retained_best_locations.csv", retained, required_retained),
        ("strain_sensor_locations.csv", sensors, required_sensor),
        ("best_feature_trials_for_every_location.csv", trials, required_trials),
    ]
    for filename, frame, required in checks:
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(
                f"{filename} is missing required columns: {sorted(missing)}"
            )

    return regions, zones, retained, trials, sensors


def zone_elements(regions: pd.DataFrame, ids) -> pd.DataFrame:
    parts = [
        p.strip()
        for p in str(ids).split(";")
        if p.strip() and p.strip().lower() != "nan"
    ]
    element_ids = [int(float(p)) for p in parts]
    frame = regions.loc[regions["elementID"].astype(int).isin(element_ids)].copy()

    if len(frame) != len(element_ids):
        found = set(frame["elementID"].astype(int))
        missing = [e for e in element_ids if e not in found]
        raise ValueError(f"Missing FE element geometry for element IDs: {missing}")

    return frame.sort_values("elementID")


def draw_segments(ax, frame: pd.DataFrame, colour, width: float):
    for row in frame.itertuples(index=False):
        ax.plot(
            [row.x1_m, row.x2_m],
            [row.y1_m, row.y2_m],
            [row.z1_m, row.z2_m],
            color=colour,
            linewidth=width,
        )


def equal_data_aspect_3d(ax):
    """
    Approximate MATLAB daspect([1 1 1]) using the displayed data ranges.
    """
    x0, x1 = ax.get_xlim3d()
    y0, y1 = ax.get_ylim3d()
    z0, z1 = ax.get_zlim3d()
    xr = max(abs(x1 - x0), 1e-12)
    yr = max(abs(y1 - y0), 1e-12)
    zr = max(abs(z1 - z0), 1e-12)
    try:
        ax.set_box_aspect((xr, yr, zr))
    except Exception:
        pass


def draw_bridge(ax, regions: pd.DataFrame, sensors: pd.DataFrame):
    """
    Python equivalent of MATLAB drawBridge().
    """
    draw_segments(ax, regions, gray(0.82), 0.35)

    xyz = sensors[["x_m", "y_m", "z_m"]].drop_duplicates().to_numpy(float)
    ax.scatter(
        xyz[:, 0],
        xyz[:, 1],
        xyz[:, 2],
        s=24,
        c=[(0.75, 0.00, 0.00)],
        edgecolors="none",
        depthshade=False,
        zorder=8,
    )

    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")
    ax.view_init(elev=21, azim=-58)
    ax.autoscale(enable=True, axis="both", tight=True)
    equal_data_aspect_3d(ax)
    ax.grid(True, alpha=0.18)


def format_probability(ax):
    ax.set_xlabel("Local stiffness loss (%)")
    ax.set_ylabel("Empirical probability (%)")
    ax.set_xlim(8, 52)
    ax.set_xticks([10, 20, 30, 40, 50])
    ax.set_ylim(0, 102)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.grid(True, alpha=0.18)


def plot_retained_locations(
    regions: pd.DataFrame,
    retained: pd.DataFrame,
    sensors: pd.DataFrame,
    figures: Path,
    dpi: int,
):
    fig = plt.figure(figsize=(13.0, 6.7))
    ax = fig.add_subplot(111, projection="3d")

    draw_bridge(ax, regions, sensors)

    for j, row in enumerate(retained.itertuples(index=False), start=1):
        q = zone_elements(regions, row.element_ids)
        draw_segments(ax, q, "black", 3.0)
        ax.text(
            row.x_m,
            row.y_m,
            row.z_m + 1.2,
            str(j),
            fontweight="bold",
            fontsize=10,
            fontfamily="DejaVu Sans",
            bbox=dict(facecolor="white", edgecolor="none", pad=1.4),
            zorder=12,
        )

    ax.set_title(
        f"Retained stiffness-loss zones ({len(retained)}); "
        "red dots: measured strain sensors"
    )

    save_figure(fig, figures / "01_retained_locations_3d", dpi)
    plt.close(fig)


def plot_retained_localisation_probabilities(
    retained: pd.DataFrame,
    trials: pd.DataFrame,
    windows: list[int],
    figures: Path,
    dpi: int,
):
    """
    Equivalent of the final/overwriting 02 plot in the MATLAB file:
    shared legend below the tiled layout.
    """
    fig, axes = plt.subplots(
        1,
        len(windows),
        figsize=(6.5 * len(windows), 5.0),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    axes = axes.ravel()

    legend_handles = []
    legend_labels = []

    # Use restrained grayscale progression so the figure remains publication-ready.
    if len(retained) > 1:
        shades = np.linspace(0.05, 0.75, len(retained))
    else:
        shades = np.array([0.15])

    markers = ["o", "s", "^", "D", "v", "P", "X", "<", ">"]

    for ii, (ax, n) in enumerate(zip(axes, windows)):
        for j, row in enumerate(retained.itertuples(index=False)):
            q = trials.loc[
                (trials["zone_index"].astype(int) == int(row.zone_index))
                & (trials["n_passages"].astype(int) == int(n))
            ].sort_values("severity_pct")

            if q.empty:
                continue

            h, = ax.plot(
                q["severity_pct"].to_numpy(float),
                100.0 * q["detected_and_correctly_localised_rate"].to_numpy(float),
                linewidth=1.3,
                marker=markers[j % len(markers)],
                markersize=4.7,
                markerfacecolor="white",
                color=gray(float(shades[j])),
                label=str(row.zone_name),
            )

            if ii == 0:
                legend_handles.append(h)
                legend_labels.append(str(row.zone_name))

        format_probability(ax)
        ax.tick_params(labelsize=16)
        ax.set_title(f"{int(n)} accumulated passages", fontsize=16)

    # Only left subplot needs y label when axes are shared.
    for ax in axes[1:]:
        ax.set_ylabel("")

    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.01),
            ncol=min(3, len(legend_handles)),
            frameon=False,
            prop={"family": "DejaVu Sans", "size": 16},
        )

    fig.subplots_adjust(bottom=0.28, wspace=0.12)
    save_figure(fig, figures / "02_retained_localisation_probabilities", dpi)
    plt.close(fig)


def plot_low_severity_trial_outcomes(
    retained: pd.DataFrame,
    trials: pd.DataFrame,
    windows: list[int],
    figures: Path,
    dpi: int,
):
    if retained.empty:
        return

    fig, axes = plt.subplots(
        1,
        len(windows),
        figsize=(6.5 * len(windows), 5.4),
        sharey=True,
        squeeze=False,
    )
    axes = axes.ravel()

    category_labels = [
        "Detected + correct",
        "Detected + wrong",
        "Correct, undetected",
        "Neither",
    ]
    shades = [0.12, 0.40, 0.65, 0.88]
    x = np.arange(len(retained))

    for ax, n in zip(axes, windows):
        q = trials.loc[
            (trials["n_passages"].astype(int) == int(n))
            & np.isclose(trials["severity_pct"].astype(float), 10.0)
        ].copy()

        q_by_id = q.set_index(q["zone_index"].astype(int), drop=False)

        ordered_rows = []
        for zone_index in retained["zone_index"].astype(int):
            if zone_index not in q_by_id.index:
                raise ValueError(
                    f"Missing 10% trial result for zone {zone_index}, {n} passages."
                )
            row = q_by_id.loc[zone_index]
            if isinstance(row, pd.DataFrame):
                if len(row) != 1:
                    raise ValueError(
                        f"Duplicate 10% trial rows for zone {zone_index}, {n} passages."
                    )
                row = row.iloc[0]
            ordered_rows.append(row)

        q = pd.DataFrame(ordered_rows).reset_index(drop=True)

        accepted = q["detected_and_correctly_localised_trials"].to_numpy(float)
        wrong = q["detected_wrong_location_trials"].to_numpy(float)
        ungated = q["correctly_localised_without_detection_trials"].to_numpy(float)
        neither = (
            q["total_trials"].to_numpy(float)
            - accepted
            - wrong
            - ungated
        )
        values = np.column_stack([accepted, wrong, ungated, neither])

        bottom = np.zeros(len(retained))
        for k in range(4):
            ax.bar(
                x,
                values[:, k],
                bottom=bottom,
                color=gray(shades[k]),
                edgecolor=gray(0.08),
                linewidth=0.45,
                width=0.72,
            )
            bottom += values[:, k]

        ax.set_xticks(x)
        ax.set_xticklabels(
            retained["zone_name"].astype(str).tolist(),
            rotation=35,
            ha="right",
        )
        ax.set_ylabel("Trial count")
        ax.grid(True, axis="y", alpha=0.18)
        ax.grid(False, axis="x")
        ax.set_title(f"10% stiffness loss; {int(n)} passages")

    handles = [
        Patch(facecolor=gray(shades[k]), edgecolor=gray(0.08), label=category_labels[k])
        for k in range(4)
    ]
    fig.legend(
        handles,
        category_labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=4,
        frameon=False,
        prop={"family": "DejaVu Sans", "size": 8},
    )

    fig.subplots_adjust(bottom=0.30, wspace=0.15)
    save_figure(fig, figures / "03_low_severity_trial_outcomes", dpi)
    plt.close(fig)


def wilson_columns(frame: pd.DataFrame, stem: str):
    lo_name = f"{stem}_wilson95_low"
    hi_name = f"{stem}_wilson95_high"
    if lo_name not in frame.columns or hi_name not in frame.columns:
        raise ValueError(
            f"Missing Wilson interval columns: {lo_name}, {hi_name}"
        )
    return (
        frame[lo_name].to_numpy(float),
        frame[hi_name].to_numpy(float),
    )


def plot_per_location_probabilities(
    zones: pd.DataFrame,
    retained: pd.DataFrame,
    trials: pd.DataFrame,
    windows: list[int],
    per_location_dir: Path,
    plot_all_locations: bool,
    dpi: int,
):
    location_ids = (
        zones["zone_index"].astype(int).tolist()
        if plot_all_locations
        else retained["zone_index"].astype(int).tolist()
    )

    metrics = [
        ("detection_rate", "detected", "--", "s", "Detection", 0.00),
        ("localisation_rate", "correctly_localised", ":", "^",
         "Correct zone (ungated)", 0.28),
        ("detected_and_correctly_localised_rate",
         "detected_and_correctly_localised", "-", "o",
         "Detected + correct zone", 0.56),
    ]

    for zone_id in location_ids:
        qz = zones.loc[zones["zone_index"].astype(int) == int(zone_id)]
        if len(qz) != 1:
            raise ValueError(f"Expected one zone-definition row for zone {zone_id}.")
        zone_name = str(qz.iloc[0]["zone_name"])

        fig, axes = plt.subplots(
            1,
            len(windows),
            figsize=(6.4 * len(windows), 4.8),
            sharex=True,
            sharey=True,
            squeeze=False,
        )
        axes = axes.ravel()

        legend_handles = []
        legend_labels = []

        for ii, (ax, n) in enumerate(zip(axes, windows)):
            q = trials.loc[
                (trials["zone_index"].astype(int) == int(zone_id))
                & (trials["n_passages"].astype(int) == int(n))
            ].sort_values("severity_pct")

            if q.empty:
                ax.set_axis_off()
                continue

            for rate_name, stem, linestyle, marker, label, shade in metrics:
                y = q[rate_name].to_numpy(float)
                low, high = wilson_columns(q, stem)
                yerr = np.vstack([
                    100.0 * np.maximum(0.0, y - low),
                    100.0 * np.maximum(0.0, high - y),
                ])

                h = ax.errorbar(
                    q["severity_pct"].to_numpy(float),
                    100.0 * y,
                    yerr=yerr,
                    linestyle=linestyle,
                    marker=marker,
                    color=gray(shade),
                    linewidth=1.2,
                    markersize=4.5,
                    capsize=2.0,
                    markerfacecolor="white",
                )

                if ii == 0:
                    legend_handles.append(h)
                    legend_labels.append(label)

            format_probability(ax)
            ax.tick_params(labelsize=14)
            ax.set_title(
                f"{int(n)} passages | {str(q.iloc[0]['feature'])}",
                fontsize=14,
            )

        fig.suptitle(
            f"{zone_name}: empirical probabilities with 95% Wilson intervals",
            fontsize=14,
            fontfamily="DejaVu Sans",
        )

        if legend_handles:
            fig.legend(
                legend_handles,
                legend_labels,
                loc="lower center",
                bbox_to_anchor=(0.5, -0.015),
                ncol=3,
                frameon=False,
                prop={"family": "DejaVu Sans", "size": 14},
            )

        fig.subplots_adjust(bottom=0.24, top=0.84, wspace=0.12)
        save_figure(
            fig,
            per_location_dir / f"zone_{int(zone_id):02d}_probabilities",
            dpi,
        )
        plt.close(fig)


# -----------------------------------------------------------------------------
# MAT-file handling for saved trial matrices.
# Supports conventional MAT files through scipy.io.loadmat, with an h5py
# fallback for MATLAB v7.3/HDF5 files.
# -----------------------------------------------------------------------------
def _load_mat_scipy(path: Path):
    try:
        from scipy.io import loadmat
    except Exception as exc:
        raise ImportError(
            "SciPy is required to read saved MATLAB trial files. "
            "Install it with: pip install scipy"
        ) from exc

    return loadmat(path, squeeze_me=True, struct_as_record=False)


def _load_mat_hdf5(path: Path):
    try:
        import h5py
    except Exception as exc:
        raise ImportError(
            "The MAT file appears to be MATLAB v7.3. Install h5py with: "
            "pip install h5py"
        ) from exc

    data = {}
    with h5py.File(path, "r") as f:
        for name in ("predictedAtom", "detected", "atomZone", "atomLoss"):
            if name not in f:
                raise KeyError(f"{name} not found in {path.name}")
            data[name] = np.array(f[name])
    return data


def load_trial_matrix(path: Path):
    try:
        data = _load_mat_scipy(path)
    except (NotImplementedError, ValueError, OSError):
        data = _load_mat_hdf5(path)

    for required in ("predictedAtom", "detected", "atomZone", "atomLoss"):
        if required not in data:
            raise KeyError(f"{required} not found in {path}")

    atom_zone = np.asarray(data["atomZone"]).reshape(-1).astype(int)
    atom_loss = np.asarray(data["atomLoss"]).reshape(-1).astype(float)

    predicted = np.asarray(data["predictedAtom"])
    detected = np.asarray(data["detected"]).astype(bool)

    predicted = np.squeeze(predicted)
    detected = np.squeeze(detected)

    if predicted.ndim == 1:
        predicted = predicted[:, None]
    if detected.ndim == 1:
        detected = detected[:, None]

    # MATLAB v7.3 arrays may appear transposed through h5py.
    n_atoms = len(atom_zone)
    if predicted.shape[0] == n_atoms and predicted.shape[1] != n_atoms:
        predicted = predicted.T
    if detected.shape[0] == n_atoms and detected.shape[1] != n_atoms:
        detected = detected.T

    if predicted.shape != detected.shape:
        raise ValueError(
            f"predictedAtom shape {predicted.shape} does not match "
            f"detected shape {detected.shape} in {path.name}"
        )

    # Stored MATLAB pipeline uses zero-based predictedAtom values.
    predicted = predicted.astype(int)

    return {
        "predictedAtom": predicted,
        "detected": detected,
        "atomZone": atom_zone,
        "atomLoss": atom_loss,
    }


def plot_spatial_predictions(
    output_root: Path,
    regions: pd.DataFrame,
    zones: pd.DataFrame,
    retained: pd.DataFrame,
    trials: pd.DataFrame,
    sensors: pd.DataFrame,
    windows: list[int],
    spatial_dir: Path,
    dpi: int,
):
    if retained.empty:
        return

    for truth in retained.itertuples(index=False):
        truth_elements = zone_elements(regions, truth.element_ids)

        for n in windows:
            rows = trials.loc[
                (trials["zone_index"].astype(int) == int(truth.zone_index))
                & (trials["n_passages"].astype(int) == int(n))
            ].sort_values("loss")

            if rows.empty:
                continue

            feature_name = str(rows.iloc[0]["feature"])
            source = output_root / "trials" / f"{feature_name}_{int(n)}pass.mat"

            if not source.is_file():
                raise FileNotFoundError(
                    f"Saved trial matrix not found:\n  {source}"
                )

            d = load_trial_matrix(source)

            fig = plt.figure(figsize=(18.0, 9.0))
            axes = [
                fig.add_subplot(2, 3, i + 1, projection="3d")
                for i in range(5)
            ]
            note_ax = fig.add_subplot(2, 3, 6)
            note_ax.axis("off")

            n_panels = min(len(rows), 5)

            for s in range(5):
                ax = axes[s]
                if s >= n_panels:
                    ax.set_axis_off()
                    continue

                row = rows.iloc[s]
                draw_bridge(ax, regions, sensors)
                draw_segments(ax, truth_elements, "black", 3.0)

                state_matches = np.where(
                    (d["atomZone"] == int(truth.zone_index))
                    & np.isclose(d["atomLoss"], float(row["loss"]), atol=1e-9, rtol=0)
                )[0]

                if len(state_matches) != 1:
                    raise ValueError(
                        "Missing or duplicate truth state in trial file "
                        f"{source.name} for zone {truth.zone_index}, loss {row['loss']}."
                    )

                a = int(state_matches[0])
                detected = d["detected"][:, a]
                predicted_atom = d["predictedAtom"][detected, a]

                if predicted_atom.size:
                    if predicted_atom.min() < 0 or predicted_atom.max() >= len(d["atomZone"]):
                        raise IndexError(
                            f"predictedAtom contains out-of-range indices in {source.name}."
                        )
                    predicted_zones = d["atomZone"][predicted_atom]
                else:
                    predicted_zones = np.array([], dtype=int)

                counts = np.array([
                    np.sum(predicted_zones == int(zone_id))
                    for zone_id in zones["zone_index"].astype(int)
                ], dtype=float)

                good = counts > 0
                n_trials = d["detected"].shape[0]
                fraction = counts / max(n_trials, 1)

                if np.any(good):
                    ax.scatter(
                        zones.loc[good, "x_m"].to_numpy(float),
                        zones.loc[good, "y_m"].to_numpy(float),
                        zones.loc[good, "z_m"].to_numpy(float),
                        s=800.0 * fraction[good],
                        facecolors="none",
                        edgecolors="black",
                        linewidths=1.3,
                        depthshade=False,
                        zorder=11,
                    )

                ax.plot(
                    [truth.x_m],
                    [truth.y_m],
                    [truth.z_m],
                    marker="x",
                    color="black",
                    markersize=12,
                    markeredgewidth=2,
                    linestyle="none",
                    zorder=12,
                )

                ax.set_title(
                    f"{int(round(float(row['severity_pct'])))}% loss | "
                    f"{int(np.sum(detected))}/{len(detected)} detected"
                )

            note_ax.text(
                0.0,
                1.0,
                "\n".join([
                    "Bold members / x: prescribed damage zone",
                    "Red dots: available measured strain sensors",
                    "Open circles: detected zone predictions",
                    "Circle area: fraction of ALL trials",
                    "Coordinates: FE zone centroids",
                    "Same-model, semi-synthetic experiment",
                ]),
                va="top",
                ha="left",
                fontfamily="DejaVu Sans",
                fontsize=11,
                transform=note_ax.transAxes,
            )

            fig.suptitle(
                f"{truth.zone_name} | {int(n)} passages | {feature_name}",
                fontsize=12,
                fontfamily="DejaVu Sans",
            )
            fig.subplots_adjust(top=0.90, wspace=0.10, hspace=0.20)

            save_figure(
                fig,
                spatial_dir
                / f"zone_{int(truth.zone_index):02d}_{int(n)}pass_predictions",
                dpi,
            )
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Redraw KW51 result figures from saved tables and trial matrices "
            "using DejaVu Sans."
        )
    )
    parser.add_argument(
        "--output-root",
        default=None,
        help=(
            "Executed KW51 output directory containing tables/ and trials/. "
            "Default: outputs/ beside this script."
        ),
    )
    parser.add_argument(
        "--retained-only",
        action="store_true",
        help=(
            "For per-location probability plots, plot only retained zones. "
            "Default is all zones, matching MATLAB plotAllLocations=true."
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=220,
        help="PNG resolution. Default: 220, matching the MATLAB script.",
    )
    args = parser.parse_args()

    output_root = resolve_output_root(args.output_root)

    figures = output_root / "figures"
    per_location = figures / "per_location"
    spatial = figures / "spatial_predictions"

    figures.mkdir(parents=True, exist_ok=True)
    per_location.mkdir(parents=True, exist_ok=True)
    spatial.mkdir(parents=True, exist_ok=True)

    regions, zones, retained, trials, sensors = read_tables(output_root)
    windows = sorted(trials["n_passages"].dropna().astype(int).unique().tolist())

    plot_retained_locations(
        regions, retained, sensors, figures, args.dpi
    )

    plot_retained_localisation_probabilities(
        retained, trials, windows, figures, args.dpi
    )

    plot_low_severity_trial_outcomes(
        retained, trials, windows, figures, args.dpi
    )

    plot_per_location_probabilities(
        zones=zones,
        retained=retained,
        trials=trials,
        windows=windows,
        per_location_dir=per_location,
        plot_all_locations=not args.retained_only,
        dpi=args.dpi,
    )

    plot_spatial_predictions(
        output_root=output_root,
        regions=regions,
        zones=zones,
        retained=retained,
        trials=trials,
        sensors=sensors,
        windows=windows,
        spatial_dir=spatial,
        dpi=args.dpi,
    )

    print(f"Saved figures under: {figures}")
    print("Font family: DejaVu Sans")


if __name__ == "__main__":
    main()
