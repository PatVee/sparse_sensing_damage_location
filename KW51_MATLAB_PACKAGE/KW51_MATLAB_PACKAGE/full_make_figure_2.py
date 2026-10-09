#!/usr/bin/env python3
r"""Rebuild KW51 figures from an ALREADY EXTRACTED output folder.

Examples
--------
    python full_make_figure.py "C:\\project\\output"
    python full_make_figure.py "C:\\project\\output" --main-only

    from full_make_figure import full_make_figure
    report = full_make_figure(r"C:\project\output", dpi=300)

Reads output/tables/*.csv and, for spatial maps, output/trials/*.mat.
Does not unzip anything, rerun a solver, choose new best features, or invent
trial outcomes. The common-family panels select feature == 'area_gram_dct8'
inside the plotting routine, in addition to location, loss and passage count.

Default outputs: output/figures/full_make_figure_python/{...}.
All original input tables and existing figures remain untouched.
Requires numpy, pandas, matplotlib and scipy. h5py is optional for v7.3 MATs.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.io import loadmat
from scipy.interpolate import PchipInterpolator

VERSION = "2026-09-24-output-folder-v1"
DEFAULT_ORDER = (
    "diagonal_A_11", "diagonal_A_14", "main_girder_A_X012",
    "main_girder_C_X084", "main_girder_C_X108",
)
KEY = ["feature", "zone_name", "n_passages", "severity_pct"]
COUNT_COLUMNS = [
    "detected_and_correctly_localised_trials", "detected_wrong_location_trials",
    "correctly_localised_without_detection_trials",
]
CATEGORIES = (
    "Detected + correct", "Detected, wrong location",
    "Correct location, not detected", "Neither",
)
CATEGORY_GRAY = (0.12, 0.46, 0.76, 1.00)
SERIES_GRAY = (0.05, 0.30, 0.48, 0.65, 0.80)
MARKERS = ("o", "s", "^", "d", "v")
LINESTYLES = ("-", "--", "-.", ":", "-")
METRICS = (
    ("detected", "detection_rate", "Detection", "--", "s", 0.10),
    ("correctly_localised", "localisation_rate", "Correct zone (ungated)", ":", "^", 0.55),
    ("detected_and_correctly_localised", "detected_and_correctly_localised_rate",
     "Detected + correct zone", "-", "o", 0.30),
)


def resolve_output_dir(output_dir: str | Path | None) -> Path:
    """Accept output, outputs, any explicitly named folder, or its tables folder."""
    if output_dir is not None:
        p = Path(output_dir).expanduser().resolve()
        if p.is_file():
            raise ValueError(f"Expected an unzipped output folder, not a file: {p}")
        if p.name.lower() == "tables" and (p / "all_location_feature_trials.csv").is_file():
            p = p.parent
        if not (p / "tables" / "all_location_feature_trials.csv").is_file():
            raise FileNotFoundError(
                f"Cannot find {p / 'tables' / 'all_location_feature_trials.csv'}. "
                "Pass the output folder itself (the folder containing tables and trials)."
            )
        return p
    here = Path(__file__).resolve().parent
    candidates: list[Path] = []
    for base in (Path.cwd(), here, here.parent):
        for candidate in (base / "output", base / "outputs", base):
            candidate = candidate.resolve()
            if candidate not in candidates:
                candidates.append(candidate)
    for p in candidates:
        if (p / "tables" / "all_location_feature_trials.csv").is_file():
            return p
    raise FileNotFoundError("No output folder found. Call full_make_figure(r'C:\\project\\output').")


def require_columns(t: pd.DataFrame, names: Iterable[str], context: str) -> None:
    missing = sorted(set(names) - set(t.columns))
    if missing:
        raise ValueError(f"{context}: missing columns: {', '.join(missing)}")


def read_csv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required result table is missing: {path}")
    t = pd.read_csv(path)
    for name in ("feature", "feature_label", "zone_name", "family", "gauge_name", "element_ids", "memberName"):
        if name in t:
            t[name] = t[name].astype(str).str.strip()
    return t


def validate_trial_table(t: pd.DataFrame, name: str) -> None:
    require_columns(t, KEY + ["zone_index", "total_trials", "detected_trials", "correctly_localised_trials"]
                    + COUNT_COLUMNS, name)
    if t.empty or t[KEY].isna().any().any():
        raise ValueError(f"{name}: empty table or missing key values.")
    if t.duplicated(KEY).any():
        sample = t.loc[t.duplicated(KEY, keep=False), KEY].head().to_dict("records")
        raise ValueError(f"{name}: duplicate feature/location/window/loss rows: {sample}")
    nums = ["total_trials", "detected_trials", "correctly_localised_trials"] + COUNT_COLUMNS
    a = t[nums].to_numpy(dtype=float)
    if not np.all(np.isfinite(a)) or np.any(a < 0) or not np.allclose(a, np.round(a), atol=1e-8, rtol=0):
        raise ValueError(f"{name}: trial counts must be finite, non-negative integers.")
    total = t.total_trials.to_numpy(dtype=float)
    counts = t[COUNT_COLUMNS].to_numpy(dtype=float)
    if np.any(total <= 0) or np.any(counts.sum(axis=1) > total):
        raise ValueError(f"{name}: invalid denominator or mutually exclusive categories exceed total_trials.")
    if not np.array_equal(t.detected_trials.to_numpy(), counts[:, 0] + counts[:, 1]):
        raise ValueError(f"{name}: detected != joint + detected_wrong.")
    if not np.array_equal(t.correctly_localised_trials.to_numpy(), counts[:, 0] + counts[:, 2]):
        raise ValueError(f"{name}: correct != joint + correct_without_detection.")
    for stem, rate, *_ in METRICS:
        require_columns(t, [rate, f"{stem}_wilson95_low", f"{stem}_wilson95_high"], name)
        if not np.allclose(t[rate], t[f"{stem}_trials"] / total, atol=1e-8, rtol=0):
            raise ValueError(f"{name}: rates do not match saved counts for {stem}.")
        lo, hi = t[f"{stem}_wilson95_low"].to_numpy(), t[f"{stem}_wilson95_high"].to_numpy()
        if not np.all(np.isfinite(lo)) or not np.all(np.isfinite(hi)) or np.any(lo < -1e-10) or np.any(hi > 1+1e-10) or np.any(lo > hi):
            raise ValueError(f"{name}: invalid interval bounds for {stem}.")
    if "loss" in t and not np.allclose(100*t.loss, t.severity_pct, atol=1e-8, rtol=0):
        raise ValueError(f"{name}: loss and severity_pct are inconsistent.")


def match_row(t: pd.DataFrame, feature: str, zone: str, passages: int, loss_pct: float) -> pd.Series:
    q = t.loc[(t.feature == feature) & (t.zone_name == zone)
              & (t.n_passages == passages) & np.isclose(t.severity_pct, loss_pct, atol=1e-8, rtol=0)]
    if len(q) != 1:
        raise ValueError(f"Expected one result for feature={feature}, zone={zone}, passages={passages}, "
                         f"loss={loss_pct:g}%; found {len(q)}.")
    return q.iloc[0]


def feature_for(t: pd.DataFrame, zone: str, passages: int) -> str:
    names = t.loc[(t.zone_name == zone) & (t.n_passages == passages), "feature"].unique()
    if len(names) != 1:
        raise ValueError(f"Expected one saved selected family for {zone}, {passages} passages; found {names}.")
    return str(names[0])


def feature_label(data: dict[str, Any], feature: str) -> str:
    q = data["catalogue"].loc[data["catalogue"].feature == feature]
    if len(q) != 1:
        raise ValueError(f"Feature catalogue must contain exactly one entry for {feature}.")
    return str(q.iloc[0].feature_label)


def zone_label(row: pd.Series) -> str:
    name = str(row.zone_name)
    m = re.fullmatch(r"diagonal_([AC])_(\d+)", name)
    if m:
        return f"{m.group(1)}{int(m.group(2)):02d} diagonal"
    m = re.fullmatch(r"main_girder_([AC])_X\d+", name)
    if m:
        return f"Main girder {m.group(1)}, x={row.x_m:.1f} m"
    return name.replace("_", " ")


def load_inputs(root: Path, common_feature: str, loss_pct: float,
                per_location: bool, location_order: Iterable[str] | None) -> dict[str, Any]:
    filenames = {
        "regions": "candidate_regions.csv", "zones": "zone_definitions.csv",
        "sensors": "strain_sensor_locations.csv", "retained": "retained_best_locations.csv",
        "selected": "retained_location_trials.csv", "all": "all_location_feature_trials.csv",
        "catalogue": "feature_catalogue.csv",
    }
    if per_location:
        filenames["best"] = "best_feature_trials_for_every_location.csv"
    data: dict[str, Any] = {key: read_csv(root / "tables" / f) for key, f in filenames.items()}
    data["source_files"] = [root / "tables" / f for f in filenames.values()]
    for key in ("all", "selected") + (("best",) if per_location else ()):
        validate_trial_table(data[key], filenames[key])
    require_columns(data["regions"], ["elementID", "family", "memberName", "x1_m", "y1_m", "z1_m", "x2_m", "y2_m", "z2_m"], "candidate_regions")
    require_columns(data["zones"], ["zone_index", "zone_name", "x_m", "y_m", "z_m", "element_ids", "physical_length_m", "n_elements"], "zone_definitions")
    require_columns(data["sensors"], ["gauge_name", "x_m", "y_m", "z_m"], "strain_sensor_locations")
    for key, col in (("zones", "zone_name"), ("zones", "zone_index"), ("regions", "elementID"), ("sensors", "gauge_name"), ("retained", "zone_name")):
        if data[key][col].duplicated().any():
            raise ValueError(f"{key}: duplicate {col}.")
    if len(data["retained"]) == 0:
        raise ValueError("No retained zones are present; the retained-zone figures cannot be made.")
    saved_names = list(data["retained"].zone_name)
    if location_order is None:
        # Preserve the five manuscript IDs; append any additional retained zones.
        order = [z for z in DEFAULT_ORDER if z in saved_names]
        order += [z for z in saved_names if z not in order]
    else:
        order = list(location_order)
        if not order or len(set(order)) != len(order) or not set(order).issubset(saved_names):
            raise ValueError("location_order must contain unique names from retained_best_locations.csv.")
    zones = data["zones"].set_index("zone_name", drop=False)
    if not set(order).issubset(zones.index):
        raise ValueError("A retained location is absent from zone_definitions.csv.")
    windows = sorted(int(x) for x in data["selected"].n_passages.unique())
    if windows != [16, 32]:
        raise ValueError(f"These manuscript layouts require passage windows [16, 32], not {windows}.")
    severities = sorted(float(x) for x in data["all"].severity_pct.unique())
    if not any(abs(s-loss_pct) < 1e-8 for s in severities):
        raise ValueError(f"Loss {loss_pct:g}% is not saved. Available: {severities}.")
    if common_feature not in set(data["catalogue"].feature):
        raise ValueError(f"Unknown common family '{common_feature}'.")
    data.update(root=root, order=order, labels={z: zone_label(zones.loc[z]) for z in order},
                windows=windows, severities=severities, common_feature=common_feature,
                loss_pct=loss_pct, zone_lookup=zones)
    thresholds = data["retained"].get("retention_threshold", pd.Series([0.5])).unique()
    if len(thresholds) != 1:
        raise ValueError("Expected one retention_threshold in retained_best_locations.csv.")
    data["retention_threshold"] = float(thresholds[0])
    # Never recompute feature rankings. Verify each saved choice against the full table.
    for key in ("selected",) + (("best",) if per_location else ()):
        for (zone, n), frame in data[key].groupby(["zone_name", "n_passages"], sort=False):
            family = feature_for(data[key], str(zone), int(n))
            if sorted(frame.severity_pct.tolist()) != severities:
                raise ValueError(f"{key}: missing/duplicate severities for {zone}, {n} passages.")
            if key == "selected":
                r = data["retained"].loc[data["retained"].zone_name == zone]
                col = f"best_feature_{int(n)}pass"
                if len(r) != 1 or col not in r or str(r.iloc[0][col]) != family:
                    raise ValueError(f"Retained feature summary disagrees with trials for {zone}, {n} passages.")
            for _, row in frame.iterrows():
                original = match_row(data["all"], family, str(zone), int(n), float(row.severity_pct))
                for col in ["total_trials", "detected_trials", "correctly_localised_trials"] + COUNT_COLUMNS:
                    if int(row[col]) != int(original[col]):
                        raise ValueError(f"{key} and full results disagree: {zone}/{n}/{row.severity_pct}/{col}.")
    for z in order:
        for n in windows:
            for loss in severities:
                match_row(data["all"], common_feature, z, n, loss)
        zone_geometry(data, z)  # Check all requested element IDs before drawing.
    return data


def zone_geometry(data: dict[str, Any], name: str) -> pd.DataFrame:
    row = data["zone_lookup"].loc[name]
    ids = [int(float(x)) for x in str(row.element_ids).split(";") if x.strip()]
    q = data["regions"].loc[data["regions"].elementID.isin(ids)].sort_values("elementID")
    if len(q) != len(ids):
        raise ValueError(f"Missing or repeated geometry elements for {name}.")
    return q


def setup_style() -> dict[str, Any]:
    available = {f.name for f in font_manager.fontManager.ttflist}
    font = "Times New Roman" if "Times New Roman" in available else "DejaVu Serif"
    return {"font.family": font, "font.size": 8.5, "axes.titlesize": 10.5,
            "axes.labelsize": 9, "axes.linewidth": 0.8, "legend.fontsize": 8.2,
            "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white"}


def style_axes(ax: Any, horizontal: bool = False) -> None:
    ax.set_facecolor("white")
    ax.set_axisbelow(True)
    ax.grid(True, axis="x" if horizontal else "both", color="0.82", alpha=0.72, linewidth=0.65)
    ax.tick_params(labelsize=8.5)


def save_figure(fig: Any, directory: Path, stem: str, dpi: int, outputs: list[str]) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    try:
        for suffix in ("pdf", "png"):
            dest = directory / f"{stem}.{suffix}"
            temporary = directory / f"{stem}.writing.{suffix}"
            fig.savefig(temporary, format=suffix, dpi=dpi, facecolor="white")
            temporary.replace(dest)
            outputs.append(str(dest))
    finally:
        plt.close(fig)


def write_plot_data(t: pd.DataFrame, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(destination, index=False, float_format="%.12g")


def draw_segments(ax: Any, frame: pd.DataFrame, gray: float, width: float) -> None:
    if frame.empty:
        return
    seg = np.stack([frame[["x1_m", "y1_m", "z1_m"]].to_numpy(),
                    frame[["x2_m", "y2_m", "z2_m"]].to_numpy()], axis=1)
    ax.add_collection3d(Line3DCollection(seg, colors=str(gray), linewidths=width))


def draw_bridge(ax: Any, data: dict[str, Any], light: bool = False) -> None:
    for family, gray, width in zip(
        ("deck_longitudinal", "deck_transverse", "main_girder", "diagonal"),
        (0.88, 0.84, 0.58, 0.67), (0.34, 0.40, 0.72, 0.78),
    ):
        frame = data["regions"].loc[data["regions"].family == family]
        draw_segments(ax, frame, 0.85 if light else gray, 0.30 if light else width)


def draw_sensors(ax: Any, data: dict[str, Any], small: bool = False) -> None:
    xyz, counts = np.unique(data["sensors"][["x_m", "y_m", "z_m"]].to_numpy(), axis=0, return_counts=True)
    ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], s=(12 if small else 34) + 10*(counts-1),
               c=[(0.75, 0.0, 0.0)], edgecolors="white", linewidths=0.65, depthshade=False)


def draw_arch_guides(ax: Any, data: dict[str, Any]) -> None:
    # Visual guide only: the CSV does not contain full arch beam connectivity.
    d = data["regions"].loc[data["regions"].family == "diagonal"]
    pts = []
    for _, member in d.groupby("memberName", sort=False):
        ends = np.unique(np.round(np.vstack([member[["x1_m", "y1_m", "z1_m"]],
                                            member[["x2_m", "y2_m", "z2_m"]]]), 8), axis=0)
        pts.append(ends[np.argmax(ends[:, 2])])
    pts = np.asarray(pts)
    for side in np.unique(np.round(pts[:, 1], 3)):
        p = pts[np.abs(pts[:, 1]-side) < 1e-3]
        p = p[np.argsort(p[:, 0])]
        p = p[np.unique(p[:, 0], return_index=True)[1]]
        if len(p) < 2:
            continue
        xf = np.linspace(p[:, 0].min(), p[:, 0].max(), 250)
        zf = PchipInterpolator(p[:, 0], p[:, 2])(xf) if len(p) >= 3 else np.interp(xf, p[:, 0], p[:, 2])
        ax.plot(xf, np.interp(xf, p[:, 0], p[:, 1]), zf, color="0.38", linewidth=1.2)


def set_bridge_axes(ax: Any, small: bool = False, view: tuple[float, float] = (20, -62)) -> None:
    ax.set(xlim=(-2, 118.5), ylim=(-1, 12.8), zlim=(0, 20),
           xticks=np.arange(0, 121, 20), yticks=[0, 6, 12], zticks=[0, 5, 10, 15, 20])
    ax.set_xlabel("x (m)" if small else "Longitudinal coordinate, x (m)", labelpad=6)
    ax.set_ylabel("y (m)", labelpad=3)
    ax.set_zlabel("z (m)" if small else "Elevation, z (m)", labelpad=2)
    ax.view_init(elev=view[0], azim=view[1])
    ax.set_box_aspect((5.4, 1.25, 1.55), zoom=1.12 if small else 1.6)
    ax.tick_params(labelsize=7 if small else 8.5, pad=0)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.fill = False
        a._axinfo["grid"]["color"] = (0.82, 0.82, 0.82, 0.60)
        a._axinfo["grid"]["linewidth"] = 0.5
    # Axes3D uses a square drawing patch even in a wide figure. With the
    # wide-bridge zoom, do not let that patch crop physically valid members.
    for artist in [*ax.lines, *ax.collections, *ax.texts]:
        artist.set_clip_on(False)


def plot_overview(data: dict[str, Any], fig_dir: Path, dpi: int, outputs: list[str], arches: bool = False) -> None:
    fig = plt.figure(figsize=(12.0, 5.2))
    ax = fig.add_axes([0.025, 0.19, 0.95, 0.70], projection="3d")
    draw_bridge(ax, data)
    if arches:
        draw_arch_guides(ax, data)
    offsets = np.array([[0, 0, 1.3], [1, 0, 2.1], [0, 0, 2.2], [-4.5, 0, 3], [0, 0, 2.2]])
    handles, labels = [], []
    for i, name in enumerate(data["order"]):
        zone = zone_geometry(data, name)
        draw_segments(ax, zone, 0.02, 4.4)
        ends = np.vstack([zone[["x1_m", "y1_m", "z1_m"]], zone[["x2_m", "y2_m", "z2_m"]]])
        ax.plot(ends[:, 0], ends[:, 1], ends[:, 2], linestyle="none", marker="o", markersize=3.7,
                markerfacecolor="white", markeredgecolor="0.02", markeredgewidth=0.65)
        centre = ends.mean(axis=0)
        tag = centre + (offsets[i] if i < len(offsets) else [0, 0, 2.2])
        tag[2] = min(tag[2], 19.3)
        ax.plot([centre[0], tag[0]], [centre[1], tag[1]], [centre[2], tag[2]-0.25], color="0.15", linewidth=0.7)
        ax.text(*tag, str(i+1), ha="center", va="center", fontsize=8.5, weight="bold",
                bbox=dict(boxstyle="circle,pad=0.30", facecolor="white", edgecolor="0.05", linewidth=0.8))
        handles.append(Line2D([], [], color="0.02", linewidth=3.2))
        labels.append(f"{i+1}  {data['labels'][name]}")
    draw_sensors(ax, data)
    handles.append(Line2D([], [], linestyle="none", marker="o", markersize=6.5,
                          markerfacecolor=(0.75, 0, 0), markeredgecolor="white"))
    labels.append("Measured strain-gauge locations")
    set_bridge_axes(ax, view=(18, -45) if arches else (20, -62))
    fig.suptitle("Retained stiffness-loss zones and strain-gauge locations", y=0.95, fontsize=11)
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.035), ncol=3, frameon=False)
    if arches:
        fig.text(0.5, 0.012, "Arch curves are visual guides through diagonal-to-arch connection coordinates, not additional FE elements.",
                 ha="center", fontsize=7)
    stem = "01b_retained_locations_3d_arch_guides" if arches else "01_retained_locations_3d"
    save_figure(fig, fig_dir, stem, dpi, outputs)


def curve_rows(data: dict[str, Any], mode: str, zone: str, passages: int) -> pd.DataFrame:
    family = data["common_feature"] if mode == "common" else feature_for(data["selected"], zone, passages)
    source = data["all"]
    q = source.loc[(source.feature == family) & (source.zone_name == zone) & (source.n_passages == passages)].copy()
    return q.sort_values("severity_pct")


def plot_joint_curves(data: dict[str, Any], fig_dir: Path, dpi: int, outputs: list[str], mode: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.55), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.072, right=0.98, bottom=0.25, top=0.81, wspace=0.09)
    handles, all_rows = [], []
    for p, n in enumerate(data["windows"]):
        ax = axes[p]
        for i, name in enumerate(data["order"]):
            q = curve_rows(data, mode, name, n)
            y = q.detected_and_correctly_localised_rate.to_numpy()*100
            lo = q.detected_and_correctly_localised_wilson95_low.to_numpy()*100
            hi = q.detected_and_correctly_localised_wilson95_high.to_numpy()*100
            h = ax.errorbar(q.severity_pct, y, yerr=np.maximum(np.vstack([y-lo, hi-y]), 0),
                            color=str(SERIES_GRAY[i % 5]), marker=MARKERS[i % 5], linestyle=LINESTYLES[i % 5],
                            markersize=4.7, markerfacecolor="white", linewidth=1.45, capsize=3)
            if p == 0:
                handles.append(h)
            q = q.copy(); q["display_location_number"] = i+1; all_rows.append(q)
        threshold = 100*data["retention_threshold"]
        ax.axhline(threshold, color="0.62", linestyle=":", linewidth=0.8)
        ax.text(max(data["severities"])-0.5, threshold+1.6, "retention threshold" if mode == "selected" else "50% reference" if threshold == 50 else f"{threshold:g}% reference",
                ha="right", va="bottom", fontsize=8, color="0.42")
        ax.set(title=f"{n} accumulated passages", xlabel="Local stiffness loss (%)",
               ylim=(0, 103), xticks=data["severities"], yticks=np.arange(0, 101, 20))
        style_axes(ax)
    axes[0].set_ylabel("Detected and correctly localised trials (%)")
    title = "Selected feature family by location and passage window" if mode == "selected" else f"Common {feature_label(data, data['common_feature'])} feature family"
    fig.suptitle(title, y=0.94, fontsize=11)
    fig.legend(handles, [data["labels"][z] for z in data["order"]], loc="lower center",
               bbox_to_anchor=(0.5, 0.025), ncol=3, frameon=False)
    stem = "02_accepted_localisation_probability" if mode == "selected" else f"05_common_{data['common_feature']}_localisation_probability"
    write_plot_data(pd.concat(all_rows, ignore_index=True), fig_dir/"plot_data"/f"{stem}.csv")
    save_figure(fig, fig_dir, stem, dpi, outputs)


def plot_low_severity_trials(data: dict[str, Any], fig_dir: Path, dpi: int,
                             outputs: list[str], mode: str) -> pd.DataFrame:
    """Every bar is explicitly matched on FEATURE + ZONE + LOSS + WINDOW."""
    all_trials = data["all"]
    display_order = list(reversed(data["order"]))
    selected_rows: list[dict[str, Any]] = []
    per_window: list[np.ndarray] = []
    for passages in data["windows"]:
        values = []
        for name in display_order:
            feature = data["common_feature"] if mode == "common" else feature_for(data["selected"], name, passages)
            # IMPORTANT: keep the explicit feature condition here; do not use
            # retained_location_trials.csv as common-family numerical input.
            row = all_trials.loc[
                (all_trials.feature == feature)
                & (all_trials.zone_name == name)
                & (all_trials.severity_pct == data["loss_pct"])
                & (all_trials.n_passages == passages)
            ]
            if len(row) != 1:
                raise ValueError(f"Expected exactly one bar source: {feature}, {name}, "
                                 f"{data['loss_pct']:g}% loss, {passages} passages; found {len(row)}.")
            r = row.iloc[0]
            three = r[COUNT_COLUMNS].to_numpy(dtype=int)
            four = np.r_[three, int(r.total_trials)-three.sum()]
            values.append(four)
            record = r.to_dict()
            record["neither_trials"] = int(four[3])
            record["display_location_number"] = data["order"].index(name)+1
            record["plot_mode"] = mode
            record["source_table"] = "tables/all_location_feature_trials.csv"
            selected_rows.append(record)
        per_window.append(np.array(values, dtype=int))
    plot_data = pd.DataFrame(selected_rows)
    totals = plot_data.total_trials.unique()
    if len(totals) != 1:
        raise ValueError("The stacked comparison requires the same trial denominator for all bars.")
    total = int(totals[0])
    if mode == "common" and set(plot_data.feature) != {data["common_feature"]}:
        raise AssertionError("The common-family plot contains more than one family.")
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.1), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.235, right=0.985, bottom=0.205, top=0.815, wspace=0.07)
    y = np.arange(1, len(display_order)+1)
    for p, passages in enumerate(data["windows"]):
        ax = axes[p]
        values = per_window[p]
        left = np.zeros(len(display_order))
        for k in range(4):
            ax.barh(y, values[:, k], left=left, height=0.64, color=str(CATEGORY_GRAY[k]),
                    edgecolor="0.08", linewidth=0.55, label=CATEGORIES[k])
            for i, count in enumerate(values[:, k]):
                if count >= 18:
                    ax.text(left[i]+count/2, y[i], str(count), ha="center", va="center", fontsize=8,
                            color="white" if k == 0 else "0.05")
            left += values[:, k]
        tick_step = 50 if total == 300 else max(1, int(math.ceil(total/6)))
        ax.set(title=f"{passages} accumulated passages", xlabel=f"Trials (out of {total})",
               xlim=(0, total), xticks=np.arange(0, total+1, tick_step), yticks=y,
               ylim=(0.5, len(display_order)+0.5))
        style_axes(ax, horizontal=True)
        if p == 0:
            ax.set_yticklabels([f"{data['order'].index(z)+1}  {data['labels'][z]}" for z in display_order])
        else:
            ax.tick_params(labelleft=False)
    loss_text = f"{data['loss_pct']:g}"
    if mode == "common":
        title = f"Prediction using common {feature_label(data, data['common_feature'])} feature family ({loss_text}% loss)"
        stem = f"04_common_{data['common_feature']}_{loss_text}pct_trial_outcomes"
    else:
        title = f"Trial outcomes at {loss_text}% local stiffness loss"
        stem = "03_low_severity_trial_outcomes"
    fig.suptitle(title, y=0.94, fontsize=11)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", bbox_to_anchor=(0.5, 0.045), ncol=4, frameon=False)
    write_plot_data(plot_data, fig_dir/"plot_data"/f"{stem}.csv")
    save_figure(fig, fig_dir, stem, dpi, outputs)
    return plot_data


def plot_per_location(data: dict[str, Any], fig_dir: Path, dpi: int, outputs: list[str]) -> None:
    directory = fig_dir / "per_location"
    for k, (_, zone) in enumerate(data["zones"].sort_values("zone_index").iterrows(), 1):
        fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.55), sharex=True, sharey=True)
        fig.subplots_adjust(left=0.065, right=0.98, bottom=0.21, top=0.78, wspace=0.09)
        source_rows, handles = [], []
        for p, n in enumerate(data["windows"]):
            ax = axes[p]
            family = feature_for(data["best"], str(zone.zone_name), n)
            q = data["all"].loc[(data["all"].feature == family) & (data["all"].zone_name == zone.zone_name)
                                & (data["all"].n_passages == n)].sort_values("severity_pct")
            source_rows.append(q)
            for stem, rate, label, ls, marker, gray in METRICS:
                y = q[rate].to_numpy()*100
                lo, hi = q[f"{stem}_wilson95_low"].to_numpy()*100, q[f"{stem}_wilson95_high"].to_numpy()*100
                h = ax.errorbar(q.severity_pct, y, yerr=np.maximum(np.vstack([y-lo, hi-y]), 0),
                                linestyle=ls, marker=marker, color=str(gray), markerfacecolor="white",
                                markersize=4, linewidth=1.25, capsize=3, label=label)
                if p == 0: handles.append(h)
            ax.set(title=f"{n} passages | {feature_label(data, family)}", xlabel="Local stiffness loss (%)",
                   ylim=(0, 103), xticks=data["severities"], yticks=np.arange(0, 101, 20))
            style_axes(ax)
        axes[0].set_ylabel("Trial fraction (%)")
        fig.suptitle(f"{zone_label(zone)}: saved selected families and 95% Wilson intervals", y=0.94, fontsize=11)
        fig.legend(handles, [m[2] for m in METRICS], loc="lower center", bbox_to_anchor=(0.5, 0.045), ncol=3, frameon=False)
        stem = f"zone_{int(zone.zone_index):02d}_probabilities"
        write_plot_data(pd.concat(source_rows, ignore_index=True), fig_dir/"plot_data"/"per_location"/f"{stem}.csv")
        save_figure(fig, directory, stem, dpi, outputs)
        if k % 10 == 0 or k == len(data["zones"]):
            print(f"  per-location plots: {k}/{len(data['zones'])}")


def load_trial_arrays(path: Path) -> dict[str, np.ndarray]:
    """Read numeric variables only; MATLAB string/MCOS objects are not needed."""
    if not path.is_file():
        raise FileNotFoundError(f"Spatial predictions need {path}. Use spatial=False or --no-spatial to omit spatial maps.")
    names = ["predictedAtom", "detected", "truthZone", "atomZone", "atomLoss", "n_passages"]
    try:
        loaded = loadmat(path, variable_names=names, squeeze_me=False)
        d = {name: np.asarray(loaded[name]) for name in names}
    except (NotImplementedError, ValueError) as exc:
        with path.open("rb") as f:
            is_v73 = b"MATLAB 7.3" in f.read(128)
        if not is_v73:
            raise ValueError(f"Cannot read numeric trial arrays from {path}: {exc}") from exc
        try:
            import h5py
        except ImportError as imp:
            raise ImportError("Install h5py to read MATLAB v7.3 trial files.") from imp
        with h5py.File(path, "r") as f:
            d = {name: np.asarray(f[name]).T for name in names}
    except KeyError as exc:
        raise ValueError(f"Missing trial variable in {path}: {exc}") from exc
    d["atomZone"] = d["atomZone"].reshape(-1).astype(np.int64)
    d["truthZone"] = d["truthZone"].reshape(-1).astype(np.int64)
    d["atomLoss"] = d["atomLoss"].reshape(-1).astype(float)
    raw_ids = d["predictedAtom"]
    if not np.all(np.isfinite(raw_ids)) or not np.array_equal(raw_ids, np.round(raw_ids)):
        raise ValueError(f"{path}: predictedAtom must contain finite integer IDs.")
    d["predictedAtom"] = raw_ids.astype(np.int64)
    if not np.all(np.isin(d["detected"], [0, 1])):
        raise ValueError(f"{path}: detected must be binary.")
    d["detected"] = d["detected"].astype(bool)
    n_atoms = len(d["atomZone"])
    if d["predictedAtom"].ndim != 2 or d["predictedAtom"].shape[1] != n_atoms or d["detected"].shape != d["predictedAtom"].shape:
        raise ValueError(f"{path}: expected trials-by-atoms matrices with {n_atoms} atom columns.")
    if len(d["atomLoss"]) != n_atoms or not np.array_equal(d["truthZone"], d["atomZone"]):
        raise ValueError(f"{path}: atom/truth maps are inconsistent.")
    # Saved predictedAtom IDs are zero-based. Do not subtract one in Python.
    if np.any(d["predictedAtom"] < 0) or np.any(d["predictedAtom"] >= n_atoms):
        raise ValueError(f"{path}: predictedAtom must contain zero-based atom IDs.")
    return d


def spatial_counts(data: dict[str, Any], d: dict[str, np.ndarray], feature: str,
                   zone_name: str, passages: int, loss_pct: float) -> tuple[pd.DataFrame, dict[str, Any]]:
    z = int(data["zone_lookup"].loc[zone_name].zone_index)
    col = np.flatnonzero((d["truthZone"] == z) & np.isclose(d["atomLoss"]*100, loss_pct, atol=1e-8, rtol=0))
    if len(col) != 1:
        raise ValueError(f"Raw MAT: expected one truth atom for zone {z}, loss {loss_pct:g}%.")
    if int(np.asarray(d["n_passages"]).item()) != passages:
        raise ValueError("MAT n_passages does not match the requested window.")
    c = int(col[0])
    predictions = d["atomZone"][d["predictedAtom"][:, c]]
    detected = d["detected"][:, c]
    correct = predictions == z
    total = len(detected)
    counts = [int(np.sum(detected & correct)), int(np.sum(detected & ~correct)),
              int(np.sum(~detected & correct)), int(np.sum(~detected & ~correct))]
    row = match_row(data["all"], feature, zone_name, passages, loss_pct)
    expected = [int(row[x]) for x in COUNT_COLUMNS]
    expected += [int(row.total_trials)-sum(expected)]
    if total != int(row.total_trials) or counts != expected:
        raise ValueError(f"MAT/CSV mismatch for {feature}/{zone_name}/{passages}/{loss_pct}: raw={counts}, CSV={expected}.")
    known = set(data["zones"].zone_index.astype(int))
    if not set(predictions).issubset(known):
        raise ValueError("A predicted zone ID has no saved zone geometry.")
    hist = data["zones"][["zone_index", "zone_name", "x_m", "y_m", "z_m"]].copy()
    hist = hist.rename(columns={"zone_index": "predicted_zone_index", "zone_name": "predicted_zone_name"})
    hist["detected_prediction_trials"] = [int(np.sum(detected & (predictions == int(k)))) for k in hist.predicted_zone_index]
    hist["fraction_of_all_trials"] = hist.detected_prediction_trials / total
    hist["truth_zone_index"] = z; hist["truth_zone_name"] = zone_name
    hist["feature"] = feature; hist["n_passages"] = passages; hist["severity_pct"] = loss_pct; hist["total_trials"] = total
    hist["is_correct_zone"] = hist.predicted_zone_index == z
    audit = dict(feature=feature, zone_name=zone_name, zone_index=z, n_passages=passages,
                 severity_pct=loss_pct, total_trials=total, detected_trials=int(detected.sum()),
                 correctly_localised_trials=int(correct.sum()),
                 detected_and_correctly_localised_trials=counts[0], detected_wrong_location_trials=counts[1],
                 correctly_localised_without_detection_trials=counts[2], neither_trials=counts[3],
                 mat_csv_counts_match=True)
    return hist, audit


def prepare_spatial(data: dict[str, Any]) -> tuple[list[dict[str, Any]], pd.DataFrame]:
    cache, prepared, audits = {}, [], []
    for mode in ("selected", "common"):
        for zone in data["order"]:
            for n in data["windows"]:
                feature = data["common_feature"] if mode == "common" else feature_for(data["selected"], zone, n)
                path = data["root"] / "trials" / f"{feature}_{n}pass.mat"
                if path not in cache:
                    cache[path] = load_trial_arrays(path)
                panels = []
                for loss in data["severities"]:
                    hist, audit = spatial_counts(data, cache[path], feature, zone, n, loss)
                    audit["plot_mode"] = mode; audit["source_mat"] = str(path)
                    panels.append((hist, audit)); audits.append(audit)
                prepared.append(dict(mode=mode, zone=zone, passages=n, feature=feature, panels=panels))
    data["source_files"] += list(cache)
    return prepared, pd.DataFrame(audits)


def plot_spatial(data: dict[str, Any], prepared: list[dict[str, Any]], fig_dir: Path, dpi: int, outputs: list[str]) -> None:
    for idx, item in enumerate(prepared, 1):
        count_panels = len(item["panels"])
        nrows = int(math.ceil((count_panels+1)/3))
        fig = plt.figure(figsize=(15.6, 4.2*nrows))
        fig.subplots_adjust(left=0.025, right=0.98, bottom=0.035, top=0.87, wspace=0.08, hspace=0.16)
        zone_name = item["zone"]
        zone = data["zone_lookup"].loc[zone_name]
        histograms = []
        for p, (hist, audit) in enumerate(item["panels"], 1):
            ax = fig.add_subplot(nrows, 3, p, projection="3d")
            draw_bridge(ax, data, light=True)
            draw_segments(ax, zone_geometry(data, zone_name), 0.02, 2.8)
            draw_sensors(ax, data, small=True)
            q = hist.loc[hist.detected_prediction_trials > 0]
            ax.scatter(q.x_m, q.y_m, q.z_m, s=850*q.fraction_of_all_trials,
                       facecolors="none", edgecolors="0.03", linewidths=1.0, depthshade=False)
            ax.scatter([zone.x_m], [zone.y_m], [zone.z_m], marker="x", s=45, color="0.0", linewidths=1.3, depthshade=False)
            set_bridge_axes(ax, small=True, view=(18, -55))
            ax.set_title(f"{audit['severity_pct']:g}% loss | {audit['detected_trials']}/{audit['total_trials']} detected",
                         fontsize=9, pad=3)
            histograms.append(hist)
        notes = fig.add_subplot(nrows, 3, count_panels+1)
        notes.axis("off")
        notes.text(0.03, 0.9,
                   "Bold member / x: prescribed stiffness-loss zone\n"
                   "Red dots: measured strain-gauge positions\n"
                   "Open circles: detected predicted zones\n"
                   "Circle area: fraction of ALL trials\n"
                   "Coordinates: saved FE zone centroids\n\n"
                   "Each MAT histogram is checked against\n"
                   "the same feature/window/loss row in the CSV.\n"
                   "Same-model, semi-synthetic assessment.",
                   fontsize=9, va="top", linespacing=1.5)
        subtitle = "Common feature" if item["mode"] == "common" else "Saved selected feature"
        fig.suptitle(f"{data['labels'][zone_name]} | {item['passages']} accumulated passages\n"
                     f"{subtitle}: {feature_label(data, item['feature'])}", y=0.96, fontsize=12)
        stem = f"zone_{int(zone.zone_index):02d}_{item['passages']}pass_predictions"
        subdir = Path("spatial_predictions") / item["mode"]
        write_plot_data(pd.concat(histograms, ignore_index=True), fig_dir/"plot_data"/subdir/f"{stem}.csv")
        save_figure(fig, fig_dir/subdir, stem, dpi, outputs)
        print(f"  spatial plots: {idx}/{len(prepared)}")


def tex_escape(value: Any) -> str:
    text = str(value)
    for ch, repl in (("&", r"\&"), ("%", r"\%"), ("_", r"\_"), ("#", r"\#")):
        text = text.replace(ch, repl)
    return text


def write_tex_table(path: Path, caption: str, label: str, headers: list[str], rows: list[list[Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = "@{}" + "l"*len(headers) + "@{}"
    lines = [r"\begin{table}[htbp]", r"\centering", r"\small", r"\setlength{\tabcolsep}{4pt}",
             f"\\caption{{{caption}}}", f"\\label{{{label}}}", f"\\begin{{tabular}}{{{cols}}}",
             r"\toprule", " & ".join(headers) + r" \\", r"\midrule"]
    lines += [" & ".join(tex_escape(x) for x in row) + r" \\" for row in rows]
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def write_tables(data: dict[str, Any], fig_dir: Path, selected_low: pd.DataFrame, common_low: pd.DataFrame) -> None:
    table_dir = fig_dir/"tables"; table_dir.mkdir(parents=True, exist_ok=True)
    feature_rows, location_rows = [], []
    for i, z in enumerate(data["order"], 1):
        zone = data["zone_lookup"].loc[z]
        location_rows.append([i, data["labels"][z], f"({zone.x_m:.1f}, {zone.y_m:.2f}, {zone.z_m:.2f})",
                              f"{zone.physical_length_m:.1f}", int(zone.n_elements)])
        for n in data["windows"]:
            f = feature_for(data["selected"], z, n)
            feature_rows.append(dict(display_location_number=i, zone_name=z, location_label=data["labels"][z],
                                     n_passages=n, feature=f, feature_label=feature_label(data, f),
                                     source_table="tables/retained_location_trials.csv"))
    features = pd.DataFrame(feature_rows)
    write_plot_data(features, table_dir/"selected_feature_by_location_and_window.csv")
    write_tex_table(table_dir/"retained_locations.tex",
                    "Retained physical stiffness-loss zones. Coordinates are saved zone centroids.", "tab:kw51-retained-zones",
                    ["ID", "Zone", "Centroid $(x,y,z)$ (m)", "Length (m)", "Elements"], location_rows)
    write_tex_table(table_dir/"selected_features.tex",
                    "Saved selected feature families, listed separately for each passage window.", "tab:kw51-selected-features",
                    ["ID", "Zone", "Passages", "Selected feature"],
                    [[r.display_location_number, r.location_label, r.n_passages, r.feature_label] for r in features.itertuples()])
    for n in data["windows"]:
        rows = []
        for i, z in enumerate(data["order"], 1):
            f = feature_for(data["selected"], z, n)
            vals = []
            for loss in data["severities"]:
                r = match_row(data["all"], f, z, n, loss)
                vals.append(f"{int(r.detected_and_correctly_localised_trials)}/{int(r.total_trials)}")
            rows.append([i, data["labels"][z]]+vals)
        write_tex_table(table_dir/f"accepted_trials_{n}.tex",
                        f"Detected and correctly localised trials for {n} accumulated passages, using saved selected families.",
                        f"tab:kw51-accepted-{n}", ["ID", "Zone"]+[f"{s:g}\\%" for s in data["severities"]], rows)
    for mode, low, filename in (("selected", selected_low, "low_severity_metrics"),
                                 ("common", common_low, "common_feature_low_severity_metrics")):
        records, rows = [], []
        for i, z in enumerate(data["order"], 1):
            for n in data["windows"]:
                r = low.loc[(low.zone_name == z) & (low.n_passages == n)].iloc[0]
                metrics = [int(r.detected_trials), int(r.correctly_localised_trials),
                           int(r.detected_and_correctly_localised_trials)]
                rows.append([i, data["labels"][z], n]+[f"{x}/{int(r.total_trials)}" for x in metrics])
                records.append(r.to_dict())
        detail = "saved selected families" if mode == "selected" else "the common "+feature_label(data, data["common_feature"])+" family"
        write_tex_table(table_dir/f"{filename}.tex",
                        f"Detection (D), ungated correct-zone selection (L), and their joint outcome at {data['loss_pct']:g}\\% stiffness loss, using {tex_escape(detail)}.",
                        f"tab:kw51-{filename.replace('_', '-')}", ["ID", "Zone", "Passages", "D", "L", "D+L"], rows)
        write_plot_data(pd.DataFrame(records), table_dir/f"{filename}.csv")
    # Export every zone for the common low-loss comparison, not just retained zones.
    q = data["all"].loc[(data["all"].feature == data["common_feature"])
                         & (data["all"].severity_pct == data["loss_pct"])].copy()
    q["shown_in_retained_zone_figures"] = q.zone_name.isin(data["order"])
    write_plot_data(q, table_dir/"common_feature_all_zones_at_low_loss.csv")


def full_make_figure(output_dir: str | Path | None = None, dpi: int = 300, *,
                     figure_dir: str | Path | None = None, per_location: bool = True,
                     spatial: bool = True, common_feature: str = "area_gram_dct8",
                     loss_pct: float = 10, location_order: Iterable[str] | None = None) -> dict[str, Any]:
    """Make the full figure set from output/tables and output/trials.

    Parameters
    ----------
    output_dir : folder containing tables/ and trials/; not a ZIP filename.
    dpi : PNG resolution (PDFs contain vector graphics).
    figure_dir : optional destination; default is below the input's figures/.
    per_location : make 68 saved-selected-family probability plots for this run.
    spatial : make selected and common-family prediction maps from saved MATs.
    common_feature : one exact saved family ID, shared by ALL displayed zones.
    loss_pct : saved loss level for the two stacked outcome plots (default 10).
    location_order : optional subset/order of retained zone names.

    Returns a manifest dict and writes it as run_manifest.json. No figures
    are left open, and no solver, classifier or feature-ranking step is run.
    """
    if isinstance(dpi, bool) or not np.isfinite(dpi) or int(dpi) != dpi or dpi <= 0:
        raise ValueError("dpi must be a positive integer.")
    if not np.isfinite(loss_pct):
        raise ValueError("loss_pct must be finite.")
    root = resolve_output_dir(output_dir)
    destination = Path(figure_dir).expanduser().resolve() if figure_dir is not None else root/"figures"/"full_make_figure_python"
    print(f"full_make_figure {VERSION}\nInput folder: {root}\nFigure folder: {destination}")
    data = load_inputs(root, common_feature, float(loss_pct), per_location, location_order)
    prepared, audit = prepare_spatial(data) if spatial else ([], pd.DataFrame())
    print(f"Loaded {len(data['all'])} feature/window/zone/loss rows; common family = {common_feature}.")
    destination.mkdir(parents=True, exist_ok=True)
    outputs: list[str] = []
    with plt.rc_context(setup_style()):
        # Common-family plot FIRST: a missing family cannot be silently skipped.
        common_low = plot_low_severity_trials(data, destination, int(dpi), outputs, "common")
        selected_low = plot_low_severity_trials(data, destination, int(dpi), outputs, "selected")
        plot_overview(data, destination, int(dpi), outputs, arches=False)
        plot_overview(data, destination, int(dpi), outputs, arches=True)
        plot_joint_curves(data, destination, int(dpi), outputs, "selected")
        plot_joint_curves(data, destination, int(dpi), outputs, "common")
        write_tables(data, destination, selected_low, common_low)
        if per_location:
            plot_per_location(data, destination, int(dpi), outputs)
        if spatial:
            plot_spatial(data, prepared, destination, int(dpi), outputs)
            write_plot_data(audit, destination/"plot_data"/"mat_csv_count_checks.csv")
    digest_rows = []
    for path in data["source_files"]:
        hasher = hashlib.sha256()
        with path.open("rb") as f:
            for block in iter(lambda: f.read(1024 * 1024), b""):
                hasher.update(block)
        digest = hasher.hexdigest()
        digest_rows.append(dict(path=str(path), sha256=digest))
    manifest = dict(version=VERSION, generated_utc=datetime.now(timezone.utc).isoformat(),
                    input_folder=str(root), figure_folder=str(destination), dpi=int(dpi),
                    common_feature=common_feature, low_loss_pct=float(loss_pct),
                    locations=data["order"], passage_windows=data["windows"], severities_pct=data["severities"],
                    all_feature_result_rows=len(data["all"]), common_plot_rows=len(common_low),
                    raw_mat_csv_checks=len(audit), figure_count=len(outputs)//2,
                    png_and_pdf_files=outputs, source_files=digest_rows,
                    interpretation="Saved exploratory semi-synthetic outcomes. A common family does not retrospectively make feature/zone selection independently validated.",
                    arch_guide_note="01b uses visual interpolation through diagonal upper ends, not saved arch-element connectivity.")
    (destination/"run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Completed: {manifest['figure_count']} figures (PNG + PDF), plus source-row CSVs and LaTeX tables.")
    print("Common-family 10%-loss source counts:" if loss_pct == 10 else "Common-family source counts:")
    print(common_low.sort_values(["n_passages", "display_location_number"])[
        ["feature", "zone_name", "n_passages", "detected_and_correctly_localised_trials", "neither_trials"]].to_string(index=False))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output_dir", nargs="?", help="Unzipped folder containing tables/ and trials/.")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--figure-dir", type=Path)
    parser.add_argument("--main-only", action="store_true", help="Six overview/comparison figures only; do not read trial MATs.")
    parser.add_argument("--no-per-location", action="store_true")
    parser.add_argument("--no-spatial", action="store_true")
    parser.add_argument("--common-feature", default="area_gram_dct8")
    parser.add_argument("--loss-pct", type=float, default=10)
    args = parser.parse_args()
    try:
        full_make_figure(args.output_dir, args.dpi, figure_dir=args.figure_dir,
                         per_location=not (args.main_only or args.no_per_location),
                         spatial=not (args.main_only or args.no_spatial),
                         common_feature=args.common_feature, loss_pct=args.loss_pct)
    except (ValueError, FileNotFoundError, ImportError, AssertionError) as exc:
        parser.exit(2, f"full_make_figure: {exc}\n")


if __name__ == "__main__":
    main()
