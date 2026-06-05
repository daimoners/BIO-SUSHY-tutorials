"""
bulk_analysis.py

Analysis helpers for the BIO-SUSHY mini-bulk tutorial.

The core tutorial-level descriptors are:
- box volume relaxation;
- mass density relaxation;
- average density over the final part of the final NPT stage.

These are intentionally simple descriptors that communicate the transition from
single-chain dynamics to material-like bulk simulations.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np


def _clean_header(header: str) -> str:
    """Normalize OpenMM StateDataReporter column headers."""
    return header.strip().strip('"').replace("#", "").strip()


def _find_column(headers: Sequence[str], candidates: Sequence[str]) -> str | None:
    lower = {h.lower(): h for h in headers}
    for candidate in candidates:
        for key, original in lower.items():
            if candidate.lower() in key:
                return original
    return None


def _read_openmm_csv(csv_file: str | os.PathLike, stage_index: int = 0, stage_name: str | None = None) -> dict:
    """Read a StateDataReporter CSV file into numeric numpy arrays."""
    path = Path(csv_file)
    if not path.exists():
        raise FileNotFoundError(f"CSV log not found: {csv_file}")

    with path.open() as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV log has no header: {csv_file}")
        original_headers = reader.fieldnames
        clean_map = {h: _clean_header(h) for h in original_headers}
        rows = []
        for row in reader:
            clean_row = {}
            for old, new in clean_map.items():
                value = row.get(old, "")
                try:
                    clean_row[new] = float(value)
                except Exception:
                    clean_row[new] = np.nan
            rows.append(clean_row)

    if not rows:
        return {"stage_index": stage_index, "stage_name": stage_name or path.stem, "n_rows": 0}

    headers = list(rows[0].keys())
    step_col = _find_column(headers, ["Step"])
    time_col = _find_column(headers, ["Time"])
    volume_col = _find_column(headers, ["Box Volume", "Volume"])
    density_col = _find_column(headers, ["Density"])
    temp_col = _find_column(headers, ["Temperature"])

    def arr(col):
        if col is None:
            return np.full(len(rows), np.nan, dtype=float)
        return np.array([row.get(col, np.nan) for row in rows], dtype=float)

    return {
        "stage_index": stage_index,
        "stage_name": stage_name or path.stem,
        "csv_file": str(path),
        "n_rows": len(rows),
        "step": arr(step_col),
        "time_ps": arr(time_col),
        "volume_nm3": arr(volume_col),
        "density_g_cm3": arr(density_col),
        "temperature_k": arr(temp_col),
    }


def _get_stage_logs_from_wm(wm) -> list[str]:
    """Find stage logs in WorkflowManager state using a few supported layouts."""
    state = getattr(wm, "state", {}) or {}

    # Preferred: wm.update_state("bulk", {"stage_logs": [...]})
    if isinstance(state.get("bulk"), dict) and state["bulk"].get("stage_logs"):
        return list(state["bulk"]["stage_logs"])

    # Alternative: wm.update_state("paths", {"bulk_stage_logs": [...]})
    if isinstance(state.get("paths"), dict) and state["paths"].get("bulk_stage_logs"):
        return list(state["paths"]["bulk_stage_logs"])

    raise ValueError(
        "Could not find bulk stage logs in workflow state. Pass stage_logs explicitly "
        "or save them with wm.update_state('bulk', {'stage_logs': [...]})"
    )


def load_bulk_timeseries(stage_logs: Sequence[str], stage_names: Sequence[str] | None = None) -> dict:
    """Load and concatenate volume/density time series from multiple stage CSV logs."""
    if not stage_logs:
        raise ValueError("stage_logs is empty.")

    stages = []
    cumulative_time = []
    volumes = []
    densities = []
    temperatures = []
    stage_ids = []
    stage_labels = []

    time_offset = 0.0
    last_time = 0.0

    for i, log in enumerate(stage_logs):
        name = stage_names[i] if stage_names and i < len(stage_names) else Path(log).stem
        data = _read_openmm_csv(log, i, name)
        if data.get("n_rows", 0) == 0:
            continue

        t = np.asarray(data["time_ps"], dtype=float)
        if not np.isfinite(t).any():
            # Fallback if time is missing: use row index.
            t = np.arange(len(data["volume_nm3"]), dtype=float)
        t = t - np.nanmin(t) + time_offset
        if len(t) > 0:
            dt = np.nanmedian(np.diff(t)) if len(t) > 1 else 0.0
            time_offset = float(np.nanmax(t) + max(dt, 0.0))
            last_time = float(np.nanmax(t))

        n = len(t)
        cumulative_time.append(t)
        volumes.append(np.asarray(data["volume_nm3"], dtype=float))
        densities.append(np.asarray(data["density_g_cm3"], dtype=float))
        temperatures.append(np.asarray(data["temperature_k"], dtype=float))
        stage_ids.append(np.full(n, i, dtype=int))
        stage_labels.extend([name] * n)
        stages.append(data)

    if not cumulative_time:
        raise ValueError("No numeric rows were found in the provided bulk logs.")

    return {
        "time_ps": np.concatenate(cumulative_time),
        "volume_nm3": np.concatenate(volumes),
        "density_g_cm3": np.concatenate(densities),
        "temperature_k": np.concatenate(temperatures),
        "stage_index": np.concatenate(stage_ids),
        "stage_label": stage_labels,
        "stages": stages,
    }


def _tail_stats(values: np.ndarray, final_fraction: float) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan, np.nan, 0
    final_fraction = float(final_fraction)
    if not 0 < final_fraction <= 1.0:
        raise ValueError("final_fraction must be > 0 and <= 1.")
    start = int(np.floor(values.size * (1.0 - final_fraction)))
    tail = values[start:]
    return float(np.mean(tail)), float(np.std(tail)), int(tail.size)


def analyze_bulk_density(
    wm=None,
    stage_logs: Sequence[str] | None = None,
    final_fraction: float = 0.5,
    save_results: bool = True,
):
    """
    Analyze density and volume relaxation from mini-bulk NPT logs.

    Parameters
    ----------
    wm
        Optional WorkflowManager. If supplied, results are saved to results.json.
    stage_logs
        List of CSV files produced by bulk_md_engine.run_bulk_annealing.
        If omitted, the function tries to read them from wm.state.
    final_fraction
        Fraction of the final stage used for final density statistics.
    save_results
        Save descriptors to WorkflowManager if available.

    Returns
    -------
    fig
        Matplotlib figure with volume and density traces.
    """
    if stage_logs is None:
        if wm is None:
            raise ValueError("Provide stage_logs or a WorkflowManager with bulk stage logs.")
        stage_logs = _get_stage_logs_from_wm(wm)

    ts = load_bulk_timeseries(stage_logs)
    final_stage = ts["stages"][-1]
    density_mean, density_std, n_density = _tail_stats(final_stage["density_g_cm3"], final_fraction)
    volume_mean, volume_std, n_volume = _tail_stats(final_stage["volume_nm3"], final_fraction)

    if wm is not None and save_results:
        wm.save_result(
            "bulk_density",
            {
                "mean": round(density_mean, 4) if np.isfinite(density_mean) else None,
                "std": round(density_std, 4) if np.isfinite(density_std) else None,
                "units": "g/cm^3",
                "analysis_window": f"last {int(final_fraction * 100)}% of final NPT stage",
                "frames_used": int(n_density),
                "note": "Mini-bulk demonstrative estimate; not a converged production density.",
            },
        )
        wm.save_result(
            "bulk_volume",
            {
                "mean": round(volume_mean, 4) if np.isfinite(volume_mean) else None,
                "std": round(volume_std, 4) if np.isfinite(volume_std) else None,
                "units": "nm^3",
                "analysis_window": f"last {int(final_fraction * 100)}% of final NPT stage",
                "frames_used": int(n_volume),
            },
        )

    print("✅ Bulk analysis complete.")
    print(f"   Final-stage density: {density_mean:.3f} ± {density_std:.3f} g/cm^3")
    print(f"   Final-stage volume:  {volume_mean:.3f} ± {volume_std:.3f} nm^3")
    print(f"   Analysis window: last {int(final_fraction * 100)}% of final NPT stage")

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    ax[0].plot(ts["time_ps"], ts["volume_nm3"], linewidth=1.5)
    ax[0].set_title("Box volume relaxation", fontsize=14)
    ax[0].set_xlabel("Cumulative reported time (ps)")
    ax[0].set_ylabel("Volume (nm$^3$)")
    ax[0].grid(True, alpha=0.3)

    ax[1].plot(ts["time_ps"], ts["density_g_cm3"], linewidth=1.5)
    if np.isfinite(density_mean):
        ax[1].axhline(density_mean, linestyle="--", alpha=0.8,
                      label=f"Final-stage mean: {density_mean:.2f} g/cm$^3$")
        ax[1].legend(loc="best")
    ax[1].set_title("Mass density relaxation", fontsize=14)
    ax[1].set_xlabel("Cumulative reported time (ps)")
    ax[1].set_ylabel("Density (g/cm$^3$)")
    ax[1].grid(True, alpha=0.3)

    # Mark stage boundaries.
    stage_index = ts["stage_index"]
    time = ts["time_ps"]
    for i in range(1, int(np.max(stage_index)) + 1):
        boundary_positions = np.where(stage_index == i)[0]
        if boundary_positions.size > 0:
            boundary_time = time[boundary_positions[0]]
            for axis in ax:
                axis.axvline(boundary_time, linestyle=":", alpha=0.5)

    fig.suptitle("Mini-bulk NPT relaxation: demonstrative material-level descriptors", y=1.03)
    plt.tight_layout()
    return fig


def print_bulk_interpretation():
    """Print a short interpretation guide for non-specialists."""
    print("🧠 Interpretation guide")
    print("- Decreasing volume usually means that the initially dilute packed box is compacting.")
    print("- Increasing density means that the polymer chains are forming a denser bulk-like system.")
    print("- The final density is a qualitative mini-bulk descriptor in this tutorial.")
    print("- Do not compare the numerical value directly with experiment unless the bulk protocol is validated and much longer.")
