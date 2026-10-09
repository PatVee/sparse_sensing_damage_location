#!/usr/bin/env python3
"""Complete Vänersborg bridge QR/virtual-sensor localisation workflow.

This script starts directly from the two raw Vänersborg ZIP archives:
    Dataset_steel_bridge_in_Sweden_1.zip
    Dataset_steel_bridge_in_Sweden_2.zip

It performs, in one reproducible pipeline:
  1. physical extraction of all raw CSV files from the two ZIP archives;
  2. strain cleaning, bridge-opening segmentation and angle alignment;
  3. assembly of 64 aligned opening records (256 samples x 12 usable gauges);
  4. physical-sensor benchmark using Area + RMS/P2P + Gram + DCT8 features;
  5. healthy POD and pivoted-QR reduction of the physical gauge set;
  6. healthy-only ridge reconstruction of the withheld strain gauges;
  7. blind virtual-sensor validation on independent baseline records;
  8. matched comparison of:
       (a) all physical sensors,
       (b) reduced QR physical sensors,
       (c) identical QR references + reconstructed virtual sensors;
  9. event-level crack/post-crack results;
 10. repeated-window bootstrap assessment;
 11. all 20 disjoint early-baseline 3/3 split sensitivity checks;
 12. CSV tables, NPZ aligned cache, PNG/PDF figures and a text summary.

The virtual channels are DATA-DRIVEN reconstructions learned from historical
healthy measurements. They are not FE-generated sensors and they are not
additional independent measurements.

Primary protocol
----------------
Usable strain gauges:
  SG1, SG2, SG3, SG7, SG8, SG9, SG10, SG11, SG13, SG14, SG15, SG16

Auxiliary raw channels:
  ch_25 : bridge-opening angle (sign reversed in the source data)
  ch_30 : temperature

Alignment:
  each opening record is baseline corrected, then interpolated on bridge angle
  from 2 to 72 degrees (128 points) and back from 72 to 2 degrees (128 points).

Event partition:
  first 3 early records : healthy fitting
  next 3 early records  : threshold calibration + blind VS validation
  documented crack file : 2023_03_09T23_45_25.csv
  damage evaluation     : crack file + the six later records

QR:
  smallest healthy POD rank retaining >=99.9% spatial energy, then pivoted QR.

Virtual transfer:
  contemporaneous standardized ridge regression with alpha=0.01, fitted only
  on the first 3 healthy records.

Spatial feature:
  signed normalized area + centered log-RMS + centered log-P2P + normalized
  Gram + first 8 orthonormal DCT-II coefficients per gauge.

Localisation truth:
  Z05_SG9_SG10, corresponding to the inspection-confirmed crack region.

No artificial damage severity or stiffness-loss coordinate is used.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import shutil
import zipfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.fft import dct
from scipy.linalg import qr
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

GAUGES = np.array([1, 2, 3, 7, 8, 9, 10, 11, 13, 14, 15, 16], dtype=int)
STRAIN_COLS = [f"ch_{g}" for g in GAUGES]
ANGLE_COL = "ch_25"
TEMP_COL = "ch_30"
USECOLS = STRAIN_COLS + [ANGLE_COL, TEMP_COL]
CRACK_FILE = "2023_03_09T23_45_25.csv"
TRUTH_ZONE = "Z05_SG9_SG10"
FIT_IDX = np.arange(3)
CAL_IDX = np.arange(3, 6)
TARGET_FS_HZ = 200.0
ALIGN_N_PER_BRANCH = 128
QR_ENERGY = 0.999
RIDGE_ALPHA = 0.01

ZONES = {
    "Z01_SG1_SG2": [1, 2],
    "Z02_SG3": [3],
    "Z03_SG7": [7],
    "Z04_SG8": [8],
    "Z05_SG9_SG10": [9, 10],
    "Z06_SG11": [11],
    "Z07_SG13_SG14": [13, 14],
    "Z08_SG15_SG16": [15, 16],
}


def sha256(path: Path, block: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def collect_zip_members(zip_paths: list[Path]) -> dict[str, tuple[Path, str]]:
    entries: dict[str, tuple[Path, str]] = {}
    for zp in zip_paths:
        with zipfile.ZipFile(zp) as zf:
            for member in zf.namelist():
                if member.lower().endswith(".csv"):
                    base = Path(member).name
                    if base in entries:
                        raise RuntimeError(f"Duplicate raw CSV basename across archives: {base}")
                    entries[base] = (zp, member)
    return dict(sorted(entries.items()))


def clean_interp(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float).copy()
    bad = (~np.isfinite(v)) | (v <= -9e5) | (np.abs(v) > 5000)
    if bad.all():
        return np.zeros_like(v)
    if bad.any():
        good = np.flatnonzero(~bad)
        bad_idx = np.flatnonzero(bad)
        if len(good) == 1:
            v[bad] = v[good[0]]
        else:
            v[bad] = np.interp(bad_idx, good, v[good])
    return v


def first_cross(a: np.ndarray, level: float, start: int = 0, direction: str = "up") -> int | None:
    if direction == "up":
        ids = np.flatnonzero((a[start:-1] < level) & (a[start + 1 :] >= level))
    else:
        ids = np.flatnonzero((a[start:-1] > level) & (a[start + 1 :] <= level))
    return start + int(ids[0]) + 1 if len(ids) else None


def interpolate_branch(angle: np.ndarray, strain: np.ndarray, i0: int, i1: int, ascending: bool,
                       nq: int = ALIGN_N_PER_BRANCH) -> np.ndarray:
    a = np.asarray(angle[i0 : i1 + 1], float)
    y = np.asarray(strain[i0 : i1 + 1, :], float)

    target = np.linspace(2.0, 72.0, nq) if ascending else np.linspace(72.0, 2.0, nq)
    if not ascending:
        a = a[::-1]
        y = y[::-1]
        target_asc = target[::-1]
    else:
        target_asc = target

    order = np.argsort(a)
    a = a[order]
    y = y[order]

    # Average repeated rounded angle bins before interpolation.
    ar = np.round(a, 2)
    unique_angle, inv = np.unique(ar, return_inverse=True)
    counts = np.bincount(inv)
    ym = np.zeros((len(unique_angle), y.shape[1]), dtype=float)
    for k in range(y.shape[1]):
        ym[:, k] = np.bincount(inv, weights=y[:, k], minlength=len(unique_angle)) / np.maximum(counts, 1)

    out = np.column_stack(
        [np.interp(target_asc, unique_angle, ym[:, k]) for k in range(y.shape[1])]
    )
    if not ascending:
        out = out[::-1]
    return out


def align_event_from_dataframe(df: pd.DataFrame, nq: int = ALIGN_N_PER_BRANCH):
    X = np.column_stack([clean_interp(df[c].to_numpy()) for c in STRAIN_COLS])
    angle = -clean_interp(df[ANGLE_COL].to_numpy())
    temp = float(np.nanmedian(clean_interp(df[TEMP_COL].to_numpy())))

    i2 = first_cross(angle, 2.0, 0, "up")
    if i2 is None:
        raise RuntimeError("No opening 2 degree crossing")
    i72 = first_cross(angle, 72.0, i2, "up")
    if i72 is None:
        raise RuntimeError("No opening 72 degree crossing")

    high = np.flatnonzero(angle >= 72.0)
    if not len(high):
        raise RuntimeError("No high-angle samples")
    j72 = int(high[-1])
    j2 = first_cross(angle, 2.0, j72, "down")
    if j2 is None:
        low = np.flatnonzero((np.arange(len(angle)) > j72) & (angle <= 4.0))
        if not len(low):
            raise RuntimeError("No closing low-angle segment")
        j2 = int(low[0])

    # Source sampling is 200 Hz. Use 7.5 s to 0.5 s before the opening crossing.
    b0 = max(0, i2 - 1500)
    b1 = max(1, i2 - 100)
    baseline = np.nanmedian(X[b0:b1, :], axis=0)
    X = X - baseline[None, :]

    opening = interpolate_branch(angle, X, i2, i72, True, nq)
    closing = interpolate_branch(angle, X, j72, j2, False, nq)
    aligned = np.vstack([opening, closing])

    audit = {
        "raw_samples": int(len(df)),
        "temperature_C": temp,
        "opening_i2": int(i2),
        "opening_i72": int(i72),
        "closing_j72": int(j72),
        "closing_j2": int(j2),
        "baseline_start": int(b0),
        "baseline_end_exclusive": int(b1),
    }
    return aligned, temp, audit


def read_raw_event(zp: Path, member: str) -> pd.DataFrame:
    with zipfile.ZipFile(zp) as zf:
        raw = zf.read(member)
    return pd.read_csv(io.BytesIO(raw), usecols=USECOLS, engine="c")


def extract_raw_archives(zip_paths: list[Path], out: Path, force: bool = False):
    """Physically extract raw CSV members from both ZIP archives.

    This stage is deliberately explicit so the reproducibility chain starts
    from the original public archives rather than an intermediate cache.
    Files are extracted under ``00_raw_extracted/<archive-stem>/``.
    """
    root = out / "00_raw_extracted"
    if force and root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True, exist_ok=True)

    extracted = []
    manifest = []
    for zp in zip_paths:
        archive_dir = root / zp.stem
        archive_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zp) as zf:
            members = [m for m in zf.namelist() if m.lower().endswith(".csv")]
            for member in members:
                dest = archive_dir / Path(member).name
                if force or not dest.exists() or dest.stat().st_size == 0:
                    with zf.open(member) as src, dest.open("wb") as dst:
                        shutil.copyfileobj(src, dst, length=1024 * 1024)
                extracted.append((zp.name, dest))
                manifest.append({
                    "archive": zp.name,
                    "member": member,
                    "file": dest.name,
                    "extracted_path": str(dest),
                    "size_bytes": dest.stat().st_size,
                })

    extracted.sort(key=lambda x: x[1].name)
    if len(extracted) != 64:
        raise RuntimeError(f"Expected 64 extracted CSV files, found {len(extracted)}")
    names = [x[1].name for x in extracted]
    if CRACK_FILE not in names:
        raise RuntimeError(f"Documented crack file not found after extraction: {CRACK_FILE}")
    pd.DataFrame(manifest).sort_values("file").to_csv(out / "00_raw_extraction_manifest.csv", index=False)
    return extracted


def _align_one_extracted(item):
    """Worker for one extracted raw CSV; returns only the compact aligned event."""
    archive_name, path_str = item
    path = Path(path_str)
    df = pd.read_csv(path, usecols=USECOLS, engine="c")
    X, temp, audit = align_event_from_dataframe(df)
    audit.update({
        "file": path.name,
        "archive": archive_name,
        "extracted_path": str(path),
    })
    return X, temp, path.name, audit


def build_aligned_dataset(extracted_files: list[tuple[str, Path]], out: Path, workers: int = 1):
    aligned, temps, files, audits = [], [], [], []
    if len(extracted_files) != 64:
        raise RuntimeError(f"Expected 64 extracted raw opening CSVs, found {len(extracted_files)}")

    tasks = [(a, str(p)) for a, p in extracted_files]
    if workers is None or workers < 1:
        workers = 1
    if workers == 1:
        iterator = map(_align_one_extracted, tasks)
        pool = None
    else:
        pool = ProcessPoolExecutor(max_workers=workers)
        iterator = pool.map(_align_one_extracted, tasks, chunksize=1)
    try:
        for k, (Xe, temp, file_name, audit) in enumerate(iterator, 1):
            aligned.append(Xe)
            temps.append(temp)
            files.append(file_name)
            audits.append(audit)
            if k % 8 == 0 or k == len(extracted_files):
                print(f"Aligned {k}/{len(extracted_files)} raw events", flush=True)
    finally:
        if pool is not None:
            pool.shutdown(wait=True)

    X = np.stack(aligned, axis=0)
    files_a = np.asarray(files)
    temps_a = np.asarray(temps, float)

    np.savez_compressed(
        out / "00_aligned_opening_cache.npz",
        aligned=X,
        files=files_a,
        gauges=GAUGES,
        temp=temps_a,
    )
    pd.DataFrame(audits).to_csv(out / "00_alignment_audit.csv", index=False)
    pd.DataFrame({
        "raw_channel": STRAIN_COLS,
        "gauge": [f"SG{g}" for g in GAUGES],
        "role": "usable strain",
    }).to_csv(out / "00_sensor_mapping.csv", index=False)
    pd.DataFrame({
        "file": files_a,
        "record_start_UTC": pd.to_datetime([f[:-4] for f in files_a], format="%Y_%m_%dT%H_%M_%S", utc=True),
        "temperature_C": temps_a,
        "role_primary": ["fit"] * 3 + ["calibration"] * 3 + ["evaluation"] * (len(files_a) - 6),
        "documented_crack_or_post": np.arange(len(files_a)) >= list(files_a).index(CRACK_FILE),
    }).to_csv(out / "00_event_inventory.csv", index=False)
    return X, files_a, GAUGES.copy(), temps_a

def area_gram_dct8(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    eps = 1e-12
    q, m = arr.shape

    meanv = arr.mean(axis=0)
    area = meanv / max(float(np.linalg.norm(meanv)), eps)

    rms = np.sqrt(np.mean(arr * arr, axis=0))
    p2p = np.ptp(arr, axis=0)
    lr = np.log(np.maximum(rms, eps))
    lp = np.log(np.maximum(p2p, eps))
    lr -= lr.mean()
    lp -= lp.mean()

    scale = math.sqrt(max(float(np.mean(arr * arr)), eps))
    y = arr / scale
    tri = np.triu_indices(m)
    gram = ((y.T @ y) / q)[tri]

    coeff = dct(y, type=2, norm="ortho", axis=0)[:8, :] / math.sqrt(q)
    dct8 = coeff.T.reshape(-1)
    return np.concatenate([area, lr, lp, gram, dct8])


def zone_scores(xin: np.ndarray, available_gauges: np.ndarray, train_idx: np.ndarray):
    available_gauges = [int(g) for g in available_gauges]
    gi = {g: i for i, g in enumerate(available_gauges)}
    zone_names, score_cols = [], []

    for name, members in ZONES.items():
        keep = [g for g in members if g in gi]
        if not keep:
            continue
        jj = [gi[g] for g in keep]
        f = np.array([area_gram_dct8(ev[:, jj]) for ev in xin])

        sd = f[train_idx].std(axis=0, ddof=1)
        positive = sd[sd > 1e-12]
        floor = (np.median(positive) * 0.05 if len(positive) else 1.0) + 1e-8
        sd = np.maximum(sd, floor)

        a = f[train_idx] / sd
        y = f / sd
        dist = ((y[:, None, :] - a[None, :, :]) ** 2).mean(axis=2)
        best = dist.min(axis=1)

        # Nearest point on every pairwise convex segment between healthy atoms.
        for i, j in combinations(range(len(train_idx)), 2):
            v = a[j] - a[i]
            t = np.clip((y - a[i]) @ v / max(float(v @ v), 1e-30), 0.0, 1.0)
            err = ((y - (a[i] + t[:, None] * v)) ** 2).mean(axis=1)
            best = np.minimum(best, err)

        zone_names.append(name)
        score_cols.append(best)

    return zone_names, np.array(score_cols).T


def select_qr_references(X: np.ndarray, gauges: np.ndarray, fit_idx: np.ndarray, energy: float = QR_ENERGY):
    snap = X[fit_idx].reshape(-1, len(gauges))
    mu = snap.mean(axis=0)
    sigma = snap.std(axis=0, ddof=1)
    sigma = np.where(sigma > 1e-12, sigma, 1.0)
    z = (snap - mu) / sigma

    _, singular, vt = np.linalg.svd(z, full_matrices=False)
    cum_energy = np.cumsum(singular ** 2) / np.sum(singular ** 2)
    rank = int(np.searchsorted(cum_energy, energy) + 1)
    phi = vt[:rank].T
    _, _, piv = qr(phi.T, pivoting=True, mode="economic")
    ref_idx = np.asarray(piv[:rank], dtype=int)
    return ref_idx, cum_energy, singular


def select_qr_fixed_count(X: np.ndarray, fit_idx: np.ndarray, k: int) -> np.ndarray:
    snap = X[fit_idx].reshape(-1, X.shape[2])
    mu = snap.mean(axis=0)
    sigma = snap.std(axis=0, ddof=1)
    sigma = np.where(sigma > 1e-12, sigma, 1.0)
    z = (snap - mu) / sigma
    _, _, vt = np.linalg.svd(z, full_matrices=False)
    phi = vt[:k].T
    _, _, piv = qr(phi.T, pivoting=True, mode="economic")
    return np.asarray(piv[:k], dtype=int)


def fit_virtual_transfer(X: np.ndarray, ref_idx: np.ndarray, withheld_idx: np.ndarray,
                         fit_idx: np.ndarray, alpha: float = RIDGE_ALPHA):
    x_train = np.vstack([X[e][:, ref_idx] for e in fit_idx])
    y_train = np.vstack([X[e][:, withheld_idx] for e in fit_idx])

    xm = x_train.mean(axis=0)
    xs = x_train.std(axis=0, ddof=1)
    ym = y_train.mean(axis=0)
    ys = y_train.std(axis=0, ddof=1)
    xs = np.where(xs > 1e-12, xs, 1.0)
    ys = np.where(ys > 1e-12, ys, 1.0)

    xz = (x_train - xm) / xs
    yz = (y_train - ym) / ys
    B = np.linalg.solve(xz.T @ xz + alpha * np.eye(xz.shape[1]), xz.T @ yz)

    xv = np.empty_like(X)
    for e in range(len(X)):
        xr = X[e][:, ref_idx]
        yp = ((xr - xm) / xs) @ B * ys + ym
        xv[e][:, ref_idx] = xr
        xv[e][:, withheld_idx] = yp

    model = {"xm": xm, "xs": xs, "ym": ym, "ys": ys, "B": B, "alpha": alpha}
    return xv, model


def validate_virtual_sensors(X: np.ndarray, Xvirtual: np.ndarray, files: np.ndarray, gauges: np.ndarray,
                             withheld_idx: np.ndarray, cal_idx: np.ndarray) -> pd.DataFrame:
    rows = []
    for e in cal_idx:
        for j in withheld_idx:
            measured = X[e, :, j]
            predicted = Xvirtual[e, :, j]
            rmse = float(np.sqrt(np.mean((predicted - measured) ** 2)))
            rows.append({
                "event": files[e],
                "gauge": f"SG{gauges[j]}",
                "NRMSE_peak_to_peak": rmse / max(float(np.ptp(measured)), 1e-12),
                "waveform_correlation": float(np.corrcoef(measured, predicted)[0, 1]),
                "amplitude_ratio": float(np.ptp(predicted) / max(float(np.ptp(measured)), 1e-12)),
            })
    return pd.DataFrame(rows)


def evaluate_configuration(label: str, xin: np.ndarray, available_gauges: np.ndarray,
                           fit_idx: np.ndarray, cal_idx: np.ndarray, damage_idx: np.ndarray,
                           crack_idx: int, representation: str):
    names, scores = zone_scores(xin, available_gauges, fit_idx)
    threshold = float(scores[cal_idx].max())
    max_score = scores.max(axis=1)
    detected = max_score > threshold
    pred = scores.argmax(axis=1)
    correct = np.array([names[i] == TRUTH_ZONE for i in pred])

    summary = {
        "configuration": label,
        "physical_sensor_count": int(len(available_gauges)),
        "representation": representation,
        "threshold": threshold,
        "detected_of_7": int(detected[damage_idx].sum()),
        "SG9_SG10_top_zone_of_7": int(correct[damage_idx].sum()),
        "joint_detect_and_localise_of_7": int((detected[damage_idx] & correct[damage_idx]).sum()),
        "median_damage_score_over_threshold": float(np.median(max_score[damage_idx] / threshold)),
        "documented_crack_detected": bool(detected[crack_idx]),
        "documented_crack_predicted_zone": names[pred[crack_idx]],
        "pre_documentation_alerts_of_51": int(detected[6:crack_idx].sum()),
    }
    return {
        "label": label,
        "zone_names": names,
        "scores": scores,
        "threshold": threshold,
        "max_score": max_score,
        "detected": detected,
        "pred": pred,
        "correct": correct,
        "summary": summary,
    }


def empirical_quantile_higher(x: np.ndarray, p: float) -> float:
    x = np.sort(np.asarray(x, float))
    idx = int(math.ceil(p * len(x)) - 1)
    idx = min(max(idx, 0), len(x) - 1)
    return float(x[idx])


def split_sensitivity(X: np.ndarray, gauges: np.ndarray, files: np.ndarray, damage_idx: np.ndarray,
                      fixed_qr_count: int, alpha: float) -> pd.DataFrame:
    rows = []
    for tr_tuple in combinations(range(6), 3):
        tr = np.array(tr_tuple, dtype=int)
        ca = np.setdiff1d(np.arange(6), tr)

        # All physical sensors.
        C = evaluate_configuration("Physical sensors only", X, gauges, tr, ca, damage_idx,
                                   damage_idx[0], "12 measured gauges")
        rows.append({
            "configuration": "Physical sensors only",
            "fit_indices": ",".join(map(str, tr)),
            "references": "",
            "detected_of_7": int(C["detected"][damage_idx].sum()),
            "correct_zone_of_7": int(C["correct"][damage_idx].sum()),
            "joint_of_7": int((C["detected"][damage_idx] & C["correct"][damage_idx]).sum()),
        })

        # Re-select a matched five-reference set for this split.
        ri = select_qr_fixed_count(X, tr, fixed_qr_count)
        rg = gauges[ri]
        Xred = X[:, :, ri]
        C = evaluate_configuration("Reduced QR physical sensors", Xred, rg, tr, ca, damage_idx,
                                   damage_idx[0], ", ".join(f"SG{g}" for g in rg))
        rows.append({
            "configuration": "Reduced QR physical sensors",
            "fit_indices": ",".join(map(str, tr)),
            "references": ",".join(f"SG{g}" for g in rg),
            "detected_of_7": int(C["detected"][damage_idx].sum()),
            "correct_zone_of_7": int(C["correct"][damage_idx].sum()),
            "joint_of_7": int((C["detected"][damage_idx] & C["correct"][damage_idx]).sum()),
        })

        wi = np.array([j for j in range(len(gauges)) if j not in ri], dtype=int)
        Xv, _ = fit_virtual_transfer(X, ri, wi, tr, alpha)
        C = evaluate_configuration("QR references + virtual sensors", Xv, gauges, tr, ca, damage_idx,
                                   damage_idx[0], f"{len(ri)} measured refs + {len(wi)} virtual gauges")
        rows.append({
            "configuration": "QR references + virtual sensors",
            "fit_indices": ",".join(map(str, tr)),
            "references": ",".join(f"SG{g}" for g in rg),
            "detected_of_7": int(C["detected"][damage_idx].sum()),
            "correct_zone_of_7": int(C["correct"][damage_idx].sum()),
            "joint_of_7": int((C["detected"][damage_idx] & C["correct"][damage_idx]).sum()),
        })
    return pd.DataFrame(rows)


def make_figures(out: Path, comparison: pd.DataFrame, validation: pd.DataFrame, event_results: pd.DataFrame,
                 split_df: pd.DataFrame, files: np.ndarray, X: np.ndarray, Xvirtual: np.ndarray,
                 gauges: np.ndarray, ref_idx: np.ndarray, withheld_idx: np.ndarray, damage_idx: np.ndarray):
    # 00: one aligned full-array event and QR role markers.
    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    event = 0
    for j, g in enumerate(gauges):
        ax.plot(np.arange(X.shape[1]), X[event, :, j], linewidth=0.9, label=f"SG{g}")
    ax.set_xlabel("Aligned opening/closing sample")
    ax.set_ylabel("Strain (source units after baseline correction)")
    ax.set_title("Vänersborg aligned baseline strain field: first reference record")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "00_aligned_baseline_strain_field.png", dpi=300)
    fig.savefig(out / "00_aligned_baseline_strain_field.pdf")
    plt.close(fig)

    # 01: matched comparison.
    fig, ax = plt.subplots(figsize=(8.7, 5.0))
    x = np.arange(len(comparison))
    width = 0.25
    ax.bar(x - width, 100 * comparison["detected_of_7"] / 7, width=width, label="Detected")
    ax.bar(x, 100 * comparison["SG9_SG10_top_zone_of_7"] / 7, width=width, label="SG9-SG10 top zone")
    ax.bar(x + width, 100 * comparison["joint_detect_and_localise_of_7"] / 7, width=width,
           label="Detected + correct")
    ax.set_xticks(x)
    ax.set_xticklabels(["All physical", "QR physical", "QR + virtual"])
    ax.set_ylim(0, 105)
    ax.set_ylabel("Fraction of seven crack/post-crack records (%)")
    ax.set_title("Vänersborg experimental localisation comparison")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.27), ncol=3)
    fig.tight_layout()
    fig.savefig(out / "01_damage_localisation_comparison.png", dpi=300)
    fig.savefig(out / "01_damage_localisation_comparison.pdf")
    plt.close(fig)

    # 02: blind virtual validation.
    gmed = validation.groupby("gauge", as_index=False)["NRMSE_peak_to_peak"].median()
    fig, ax = plt.subplots(figsize=(8.7, 5.0))
    ax.bar(gmed["gauge"], 100 * gmed["NRMSE_peak_to_peak"])
    ax.set_ylabel("Median NRMSE (% of measured peak-to-peak)")
    ax.set_title("Blind validation of withheld virtual strain gauges")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "02_virtual_sensor_validation_NRMSE.png", dpi=300)
    fig.savefig(out / "02_virtual_sensor_validation_NRMSE.pdf")
    plt.close(fig)

    # 03: eventwise detection margin.
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    xpos = np.arange(len(damage_idx))
    for label in comparison["configuration"]:
        q = event_results[event_results["configuration"] == label]
        ax.plot(xpos, q["score_over_threshold"], marker="o", label=label)
    ax.axhline(1.0, linestyle="--", linewidth=1.1)
    ax.set_xticks(xpos)
    ax.set_xticklabels(["Crack"] + [f"Post {i}" for i in range(1, len(damage_idx))])
    ax.set_ylabel("Damage score / threshold")
    ax.set_title("Eventwise detection margin")
    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "03_eventwise_detection_margin.png", dpi=300)
    fig.savefig(out / "03_eventwise_detection_margin.pdf")
    plt.close(fig)

    # 04: baseline-split sensitivity.
    sens = split_df.groupby("configuration")["joint_of_7"].agg(["min", "median", "max"]).reset_index()
    desired = ["Physical sensors only", "Reduced QR physical sensors", "QR references + virtual sensors"]
    sens["order"] = sens["configuration"].map({x: i for i, x in enumerate(desired)})
    sens = sens.sort_values("order")
    fig, ax = plt.subplots(figsize=(8.7, 5.0))
    x = np.arange(len(sens))
    ax.errorbar(
        x,
        100 * sens["median"].to_numpy() / 7,
        yerr=np.vstack([
            100 * (sens["median"] - sens["min"]).to_numpy() / 7,
            100 * (sens["max"] - sens["median"]).to_numpy() / 7,
        ]),
        fmt="o",
        capsize=5,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(["All physical", "QR physical", "QR + virtual"])
    ax.set_ylim(-3, 105)
    ax.set_ylabel("Detected + correctly localised (%)")
    ax.set_title("Sensitivity to all 20 early-baseline 3/3 splits")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "04_baseline_split_sensitivity.png", dpi=300)
    fig.savefig(out / "04_baseline_split_sensitivity.pdf")
    plt.close(fig)

    # 05: representative real-vs-virtual waveform for the first withheld gauge on calibration record 4.
    e = int(CAL_IDX[0])
    j = int(withheld_idx[0])
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ax.plot(X[e, :, j], label=f"Measured SG{gauges[j]}")
    ax.plot(Xvirtual[e, :, j], linestyle="--", label=f"Virtual SG{gauges[j]}")
    ax.set_xlabel("Aligned opening/closing sample")
    ax.set_ylabel("Strain")
    ax.set_title("Representative blind virtual-strain validation")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "05_virtual_sensor_waveform_example.png", dpi=300)
    fig.savefig(out / "05_virtual_sensor_waveform_example.pdf")
    plt.close(fig)


def make_conceptual_figures(out: Path, comparison: pd.DataFrame, validation: pd.DataFrame,
                            split_df: pd.DataFrame, configs: list[dict], crack_idx: int):
    """Create manuscript-oriented conceptual result figures in PNG and PDF."""

    def save(fig, stem):
        fig.tight_layout()
        fig.savefig(out / f"{stem}.png", dpi=300, bbox_inches="tight")
        fig.savefig(out / f"{stem}.pdf", bbox_inches="tight")
        plt.close(fig)

    # 06: conceptual information-recovery pathway.
    fig, ax = plt.subplots(figsize=(11.5, 5.2))
    ax.set_xlim(0, 12); ax.set_ylim(0, 5.2); ax.axis("off")
    boxes = [
        (0.4, 1.55, 2.2, 2.0, "12 physical\nstrain gauges", f"{int(comparison.iloc[0].joint_detect_and_localise_of_7)}/7 joint\nlocalisation"),
        (3.3, 1.55, 2.2, 2.0, f"QR selection\n{int(comparison.iloc[1].physical_sensor_count)} references", f"{int(comparison.iloc[1].joint_detect_and_localise_of_7)}/7 joint\nafter reduction"),
        (6.2, 1.55, 2.2, 2.0, "Healthy-trained\nvirtual transfer", f"{12-int(comparison.iloc[1].physical_sensor_count)} withheld channels\nreconstructed"),
        (9.1, 1.55, 2.2, 2.0, f"{int(comparison.iloc[1].physical_sensor_count)} measured +\n{12-int(comparison.iloc[1].physical_sensor_count)} virtual", f"{int(comparison.iloc[2].joint_detect_and_localise_of_7)}/7 joint\nlocalisation"),
    ]
    for x0, y0, ww, hh, title, result in boxes:
        ax.add_patch(FancyBboxPatch((x0, y0), ww, hh, boxstyle="round,pad=0.03"))
        ax.text(x0+ww/2, y0+1.25, title, ha="center", va="center", fontsize=12, fontweight="bold")
        ax.text(x0+ww/2, y0+0.52, result, ha="center", va="center", fontsize=10)
    for a, b in [(2.6,3.3),(5.5,6.2),(8.4,9.1)]:
        ax.add_patch(FancyArrowPatch((a,2.55),(b,2.55),arrowstyle="->",mutation_scale=14))
    ax.text(6.0,0.45,
        "Virtual channels are historically calibrated reconstructions of removed physical gauges;\n"
        "they are not independent measurements and are not FE-generated sensors.",
        ha="center", va="center", fontsize=10)
    ax.set_title("Vänersborg sensing experiment: information loss and recovery", pad=14)
    save(fig, "06_information_recovery_pathway")

    # 07: spatial crack-zone fingerprint at the documented crack event.
    all_zone_names = list(ZONES.keys())
    M = np.full((3, len(all_zone_names)), np.nan, dtype=float)
    for r, C in enumerate(configs):
        for j, name in enumerate(C["zone_names"]):
            M[r, all_zone_names.index(name)] = C["scores"][crack_idx, j]
    rowmax = np.nanmax(M, axis=1, keepdims=True)
    Mnorm = M / np.maximum(rowmax, 1e-30)
    fig, ax = plt.subplots(figsize=(11.0, 4.6))
    im = ax.imshow(Mnorm, aspect="auto")
    ax.set_yticks([0,1,2]); ax.set_yticklabels(["All physical (12)","QR physical (5)","QR + virtual (5+7)"])
    short = ["SG1–2","SG3","SG7","SG8","SG9–10","SG11","SG13–14","SG15–16"]
    ax.set_xticks(np.arange(8)); ax.set_xticklabels(short)
    ax.set_xlabel("Sensor zone"); ax.set_title("Spatial anomaly fingerprint at the documented crack record")
    cb = fig.colorbar(im, ax=ax); cb.set_label("Zone score normalized within each sensing configuration")
    for i in range(3):
        for j in range(8):
            txt = "N/A" if not np.isfinite(Mnorm[i,j]) else f"{Mnorm[i,j]:.2f}"
            ax.text(j,i,txt,ha="center",va="center",fontsize=9)
    ax.add_patch(Rectangle((4-0.5,-0.5),1,3,fill=False,linewidth=2.2))
    ax.text(4,-0.82,"inspection truth",ha="center",va="bottom",fontsize=9)
    save(fig, "07_crack_zone_fingerprint")

    # 08: virtual-sensor validation trust map.
    gmed = validation.groupby("gauge", as_index=False).agg(
        median_NRMSE=("NRMSE_peak_to_peak","median"),
        median_corr=("waveform_correlation","median"),
        median_amp=("amplitude_ratio","median"),
    )
    fig, ax = plt.subplots(figsize=(7.5,5.6))
    ax.scatter(100*gmed["median_NRMSE"], gmed["median_corr"], s=65)
    for _, rr in gmed.iterrows():
        ax.annotate(rr["gauge"], (100*rr["median_NRMSE"], rr["median_corr"]),
                    xytext=(4,3), textcoords="offset points", fontsize=9)
    ax.set_xlabel("Median NRMSE (% of measured peak-to-peak strain)")
    ax.set_ylabel("Median waveform correlation")
    ax.set_title("Blind validation of reconstructed strain channels")
    ax.grid(True, alpha=0.25)
    save(fig, "08_virtual_sensor_validation_map")

    # 09: all 20 baseline 3/3 partitions, not only min/median/max.
    fig, ax = plt.subplots(figsize=(9.2,5.6))
    labels = ["Physical sensors only","Reduced QR physical sensors","QR references + virtual sensors"]
    xpos = np.arange(3)
    for i, label in enumerate(labels):
        vals = split_df.loc[split_df["configuration"]==label,"joint_of_7"].to_numpy()
        offsets = np.linspace(-0.17,0.17,len(vals))
        ax.scatter(np.full(len(vals),xpos[i])+offsets, vals, s=34)
        ax.plot([xpos[i]-0.22,xpos[i]+0.22],[np.median(vals),np.median(vals)],linewidth=2.0)
    ax.set_xticks(xpos); ax.set_xticklabels(["All physical\n(12)","QR physical\n(5)","QR + virtual\n(5+7)"])
    ax.set_ylim(-0.3,7.5); ax.set_yticks(range(8))
    ax.set_ylabel("Detected + correctly localised records (out of 7)")
    ax.set_title("Robustness across all 20 early-baseline 3/3 partitions")
    ax.grid(axis="y", alpha=0.25)
    save(fig, "09_baseline_partition_robustness")

    # Table of the documented-crack zone scores used in Figure 07.
    rows=[]
    for r, C in enumerate(configs):
        for j, z in enumerate(all_zone_names):
            rows.append({
                "configuration": C["label"],
                "zone": z,
                "raw_score": M[r,j],
                "normalized_score": Mnorm[r,j],
                "inspection_truth": z == TRUTH_ZONE,
            })
    pd.DataFrame(rows).to_csv(out / "10_documented_crack_zone_scores.csv", index=False)

    pd.DataFrame({
        "figure": [
            "01_damage_localisation_comparison",
            "06_information_recovery_pathway",
            "07_crack_zone_fingerprint",
            "08_virtual_sensor_validation_map",
            "09_baseline_partition_robustness",
        ],
        "role": [
            "Primary quantitative result",
            "Conceptual contribution",
            "Spatial localisation evidence",
            "Virtual-sensor credibility",
            "Robustness / limitation",
        ]
    }).to_csv(out / "11_figure_manifest.csv", index=False)


def run_analysis(X: np.ndarray, files: np.ndarray, gauges: np.ndarray, temp: np.ndarray, out: Path):
    crack_idx = list(files).index(CRACK_FILE)
    damage_idx = np.arange(crack_idx, len(files))

    ref_idx, cum_energy, singular = select_qr_references(X, gauges, FIT_IDX, QR_ENERGY)
    ref_gauges = gauges[ref_idx]
    withheld_idx = np.array([j for j in range(len(gauges)) if j not in ref_idx], dtype=int)
    withheld_gauges = gauges[withheld_idx]
    Xvirtual, virtual_model = fit_virtual_transfer(X, ref_idx, withheld_idx, FIT_IDX, RIDGE_ALPHA)

    validation = validate_virtual_sensors(X, Xvirtual, files, gauges, withheld_idx, CAL_IDX)
    validation.to_csv(out / "03_virtual_sensor_blind_validation.csv", index=False)

    C1 = evaluate_configuration(
        "Physical sensors only", X, gauges, FIT_IDX, CAL_IDX, damage_idx, crack_idx, "12 measured gauges"
    )
    Xred = X[:, :, ref_idx]
    C2 = evaluate_configuration(
        "Reduced QR physical sensors", Xred, ref_gauges, FIT_IDX, CAL_IDX, damage_idx, crack_idx,
        ", ".join(f"SG{g}" for g in ref_gauges),
    )
    C3 = evaluate_configuration(
        "QR references + virtual sensors", Xvirtual, gauges, FIT_IDX, CAL_IDX, damage_idx, crack_idx,
        f"{len(ref_gauges)} measured refs + {len(withheld_gauges)} virtual gauges",
    )
    C1["summary"]["physical_sensor_count"] = 12
    C2["summary"]["physical_sensor_count"] = len(ref_gauges)
    C3["summary"]["physical_sensor_count"] = len(ref_gauges)
    configs = [C1, C2, C3]

    comparison = pd.DataFrame([c["summary"] for c in configs])
    comparison.to_csv(out / "01_three_configuration_summary.csv", index=False)

    event_rows = []
    for C in configs:
        for i in damage_idx:
            event_rows.append({
                "configuration": C["label"],
                "file": files[i],
                "is_documented_crack_record": bool(i == crack_idx),
                "score": float(C["max_score"][i]),
                "threshold": float(C["threshold"]),
                "score_over_threshold": float(C["max_score"][i] / C["threshold"]),
                "detected": bool(C["detected"][i]),
                "predicted_zone": C["zone_names"][C["pred"][i]],
                "correct_SG9_SG10_zone": bool(C["correct"][i]),
                "joint_correct": bool(C["detected"][i] and C["correct"][i]),
            })
    event_results = pd.DataFrame(event_rows)
    event_results.to_csv(out / "02_crack_and_post_event_results.csv", index=False)

    qr_table = pd.DataFrame({
        "QR_rank": np.arange(1, len(ref_gauges) + 1),
        "physical_reference": [f"SG{g}" for g in ref_gauges],
        "selection_rule": "smallest POD rank reaching >=99.9% fit-record spatial energy",
        "cumulative_energy_at_selected_rank": cum_energy[len(ref_gauges) - 1],
    })
    qr_table.to_csv(out / "06_QR_reference_selection.csv", index=False)
    pd.DataFrame({
        "POD_rank": np.arange(1, len(cum_energy) + 1),
        "singular_value": singular,
        "cumulative_healthy_spatial_energy": cum_energy,
    }).to_csv(out / "07_healthy_spatial_POD_energy.csv", index=False)

    # Save frozen virtual transfer matrices/normalisation explicitly.
    np.savez_compressed(
        out / "08_frozen_virtual_transfer_model.npz",
        reference_indices=ref_idx,
        reference_gauges=ref_gauges,
        withheld_indices=withheld_idx,
        withheld_gauges=withheld_gauges,
        xm=virtual_model["xm"],
        xs=virtual_model["xs"],
        ym=virtual_model["ym"],
        ys=virtual_model["ys"],
        ridge_B=virtual_model["B"],
        ridge_alpha=RIDGE_ALPHA,
    )

    # Conditional repeated-window resampling.
    rng = np.random.default_rng(20261005)
    rep_rows = []
    for C in configs:
        for window in [1, 3, 5]:
            cal_draw = rng.choice(CAL_IDX, size=(5000, window), replace=True)
            cal_acc = C["scores"][cal_draw].sum(axis=1)
            threshold_acc = empirical_quantile_higher(cal_acc.max(axis=1), 0.95)

            damage_draw = rng.choice(damage_idx, size=(1000, window), replace=True)
            dmg_acc = C["scores"][damage_draw].sum(axis=1)
            det = dmg_acc.max(axis=1) > threshold_acc
            pp = dmg_acc.argmax(axis=1)
            cc = np.array([C["zone_names"][i] == TRUTH_ZONE for i in pp])
            rep_rows.append({
                "configuration": C["label"],
                "window_events": window,
                "trials": 1000,
                "accumulated_threshold": threshold_acc,
                "detection_rate": float(det.mean()),
                "SG9_SG10_selection_rate": float(cc.mean()),
                "joint_rate": float((det & cc).mean()),
            })
    repeatability = pd.DataFrame(rep_rows)
    repeatability.to_csv(out / "04_repeatability_bootstrap.csv", index=False)

    split_df = split_sensitivity(X, gauges, files, damage_idx, len(ref_idx), RIDGE_ALPHA)
    split_df.to_csv(out / "05_baseline_split_sensitivity.csv", index=False)

    # Summary of virtual validation and robustness.
    split_summary = split_df.groupby("configuration")["joint_of_7"].agg(["min", "median", "max"])
    key = {
        "n_events": int(len(files)),
        "aligned_samples_per_event": int(X.shape[1]),
        "usable_physical_gauges": [int(x) for x in gauges],
        "fit_files": [str(files[i]) for i in FIT_IDX],
        "calibration_files": [str(files[i]) for i in CAL_IDX],
        "documented_crack_file": CRACK_FILE,
        "damage_test_files": [str(files[i]) for i in damage_idx],
        "truth_zone": TRUTH_ZONE,
        "qr_energy_target": QR_ENERGY,
        "qr_rank": int(len(ref_idx)),
        "qr_reference_gauges": [int(x) for x in ref_gauges],
        "virtual_gauges": [int(x) for x in withheld_gauges],
        "ridge_alpha": RIDGE_ALPHA,
        "virtual_validation_median_NRMSE": float(validation["NRMSE_peak_to_peak"].median()),
        "virtual_validation_median_correlation": float(validation["waveform_correlation"].median()),
        "virtual_validation_median_amplitude_ratio": float(validation["amplitude_ratio"].median()),
        "primary_results": comparison.to_dict(orient="records"),
        "split_sensitivity_joint_counts": {
            idx: {k: float(v) for k, v in row.items()}
            for idx, row in split_summary.to_dict(orient="index").items()
        },
        "interpretation": (
            "The virtual channels are historically calibrated data-driven reconstructions from the QR-selected "
            "physical gauges. They are not FE-generated sensors and are not independent measurements."
        ),
    }
    (out / "09_run_summary.json").write_text(json.dumps(key, indent=2), encoding="utf-8")

    make_figures(out, comparison, validation, event_results, split_df, files, X, Xvirtual,
                 gauges, ref_idx, withheld_idx, damage_idx)
    make_conceptual_figures(out, comparison, validation, split_df, configs, crack_idx)

    summary_text = f"""VÄNERSBORG FULL PIPELINE EXECUTED SUMMARY
=========================================
Raw events: {len(files)}
Aligned shape: {X.shape}
Documented crack record: {CRACK_FILE}
Truth zone: {TRUTH_ZONE}

QR selection
------------
Energy target: {100*QR_ENERGY:.2f}%
Selected rank: {len(ref_gauges)}
References: {', '.join(f'SG{g}' for g in ref_gauges)}
Virtual/withheld gauges: {', '.join(f'SG{g}' for g in withheld_gauges)}

Blind virtual-sensor validation
-------------------------------
Median NRMSE: {100*validation['NRMSE_peak_to_peak'].median():.2f}% of measured peak-to-peak
Median correlation: {validation['waveform_correlation'].median():.4f}
Median amplitude ratio: {validation['amplitude_ratio'].median():.4f}

Primary three-way comparison
----------------------------
{comparison[['configuration','detected_of_7','SG9_SG10_top_zone_of_7','joint_detect_and_localise_of_7']].to_string(index=False)}

Interpretation boundary
-----------------------
The virtual gauges are data-driven reconstructions learned exclusively from historical healthy
full-array measurements. They are not FE-generated sensors and they do not add independent
information beyond the physical references. The primary experiment asks whether historical
full-array calibration can preserve localisation after permanent sensor reduction.
"""
    (out / "README_RESULTS.txt").write_text(summary_text, encoding="utf-8")
    print(summary_text)
    return comparison


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zip1", type=Path, required=True, help="Dataset_steel_bridge_in_Sweden_1.zip")
    ap.add_argument("--zip2", type=Path, required=True, help="Dataset_steel_bridge_in_Sweden_2.zip")
    ap.add_argument("--out", type=Path, default=Path("vanersborg_full_results"), help="Output directory")
    ap.add_argument("--reuse-cache", action="store_true", help="Reuse 00_aligned_opening_cache.npz if present")
    ap.add_argument("--reuse-extracted", action="store_true", help="Reuse 00_raw_extracted CSVs if already present")
    ap.add_argument("--workers", type=int, default=4, help="Parallel CSV alignment workers (default: 4)")
    args = ap.parse_args()

    for p in [args.zip1, args.zip2]:
        if not p.is_file():
            raise FileNotFoundError(p)
    args.out.mkdir(parents=True, exist_ok=True)

    input_manifest = pd.DataFrame([
        {"file": args.zip1.name, "size_bytes": args.zip1.stat().st_size, "sha256": sha256(args.zip1)},
        {"file": args.zip2.name, "size_bytes": args.zip2.stat().st_size, "sha256": sha256(args.zip2)},
    ])
    input_manifest.to_csv(args.out / "00_input_manifest.csv", index=False)

    cache = args.out / "00_aligned_opening_cache.npz"
    if args.reuse_cache and cache.exists():
        d = np.load(cache, allow_pickle=True)
        X = d["aligned"].astype(float)
        files = d["files"].astype(str)
        gauges = d["gauges"].astype(int)
        temp = d["temp"].astype(float)
        print(f"Reused aligned cache: {cache}")
    else:
        extracted = extract_raw_archives([args.zip1, args.zip2], args.out, force=not args.reuse_extracted)
        X, files, gauges, temp = build_aligned_dataset(extracted, args.out, workers=args.workers)

    if X.shape != (64, 256, 12):
        raise RuntimeError(f"Unexpected aligned array shape {X.shape}; expected (64, 256, 12)")
    run_analysis(X, files, gauges, temp, args.out)


if __name__ == "__main__":
    main()
