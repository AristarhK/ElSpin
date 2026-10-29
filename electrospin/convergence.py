"""Independent field, time-step and parcel-resolution sensitivity studies.

This reports changes in numerical results, not a convergence certificate or an
experimental validation.  All physical parameters (including perturbation and
Coulomb softening) remain fixed.  Finer parcel intervals change quadrature of a
fixed physical inlet; they do not change its flow rate or current.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np

from .config import ROOT, load_config, write_config
from .field import solve_field


FIELD_METRICS = (
    "needle_charge_C", "collector_charge_C", "chamber_charge_C",
    "tip_axis_ez_V_m", "axis_ez_at_0p01_gap_V_m", "axis_ez_at_0p1_gap_V_m",
    "axis_ez_at_0p5_gap_V_m", "relative_linear_residual",
    "relative_boundary_charge_imbalance",
)
DYNAMIC_METRICS = (
    "deposited_parcels", "first_impact_time_s", "rms_impact_radius_material_m",
    "mean_wet_diameter_at_hits_m", "mean_dry_diameter_same_length_at_hits_m",
    "delivered_polymer_mass_kg", "delivered_charge_C", "mapped_centerline_length_m",
    "orientation_nematic_order_0_to_1", "final_max_radial_excursion_m",
    "maximum_mass_relative_error", "maximum_polymer_relative_error",
    "maximum_charge_relative_error", "sampled_min_dt_s", "sampled_max_dt_s",
)


def _field_metrics(directory: Path, config: dict) -> dict:
    summary = json.loads((directory / "field_summary.json").read_text())
    data = np.genfromtxt(directory / "potential_rz.csv", delimiter=",", names=True)
    axis = data[data["r_m"] == data["r_m"].min()]
    axis.sort(order="z_m")
    samples = np.interp(np.array([0., .01, .1, .5]) * config["gap"],
                        axis["z_m"], axis["E_z_V_m"])
    return {
        **{f"{electrode}_charge_C": charge
           for electrode, charge in summary["electrode_charge_C"].items()},
        "tip_axis_ez_V_m": float(samples[0]),
        "axis_ez_at_0p01_gap_V_m": float(samples[1]),
        "axis_ez_at_0p1_gap_V_m": float(samples[2]),
        "axis_ez_at_0p5_gap_V_m": float(samples[3]),
        "relative_linear_residual": summary["relative_linear_residual"],
        "relative_boundary_charge_imbalance": summary["relative_boundary_charge_imbalance"],
    }


def _dynamic_metrics(summary: dict, directory: Path) -> dict:
    result = {name: summary.get(name) for name in DYNAMIC_METRICS}
    data = np.genfromtxt(directory / "diagnostics.csv", delimiter=",", names=True)
    if np.size(data):
        values = np.atleast_1d(data["dt_s"])
        values = values[values > 0]
        if len(values):
            result["sampled_min_dt_s"] = float(values.min())
            result["sampled_max_dt_s"] = float(values.max())
    return result


def _save_table(out: Path, rows: list[dict]) -> None:
    # A planned but unexecuted/failed case has no numerical estimates. Never
    # substitute zero, reuse another result, or silently omit the case.
    names = list(dict.fromkeys(name for row in rows for name in row))
    with (out / "convergence.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, names)
        writer.writeheader()
        writer.writerows(rows)


def _changes(current: dict, previous: dict | None, metrics: tuple[str, ...]) -> None:
    if previous is None or previous["status"] != "COMPLETED":
        return
    for metric in metrics:
        a, b = current.get(metric), previous.get(metric)
        if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
            continue
        if not math.isfinite(a) or not math.isfinite(b):
            continue
        current[metric + "_absolute_change_previous"] = abs(a - b)
        # No relative estimate for an exactly zero current value unless both
        # values vanish. Physical zeros must not produce infinity in JSON/CSV.
        current[metric + "_relative_change_previous"] = (abs(a-b)/abs(a) if a else (0. if not b else None))


def study(config: dict, out: Path, *, field_only: bool = False, debug: bool = False) -> int:
    out = Path(out).resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Output directory is not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    write_config(config, out / "base_config.json")
    designs = []
    for level, nr, nz in (("coarse", 48, 80), ("standard", 96, 160), ("fine", 192, 320)):
        designs.append(("field", level, {**config, "field_nr": nr, "field_nz": nz}))
    for kind, key in (("time_step", "dt_max"), ("parcel_interval", "injection_interval")):
        for level, factor in (("coarse", 2.), ("standard", 1.), ("fine", .5)):
            case = {**config, key: config[key] * factor}
            case["max_beads"] = max(config["max_beads"], math.ceil(case["duration"] / case["injection_interval"]) + 4)
            designs.append((kind, level, case))
    rows = [dict(study=kind, level=level, status="UNRUN", directory=str(out / kind / level),
                 field_nr=case["field_nr"], field_nz=case["field_nz"],
                 dt_max=case["dt_max"], injection_interval=case["injection_interval"],
                 reason="excluded_by_field_only" if field_only and kind != "field" else "pending")
            for kind, level, case in designs]
    _save_table(out, rows)
    instructions = {
        "purpose": "Numerical sensitivity; no automatic convergence or physical-validation claim.",
        "field": "48x80, 96x160, 192x320 meshes of identical needle/meniscus geometry; local tip field has staircase sensitivity.",
        "time_step": "dt_max factors 2, 1, 0.5; adaptive stability limits can make actual steps smaller than all three maxima.",
        "parcel_interval": "injection_interval factors 2, 1, 0.5 at constant flow/current; capacity max_beads increases when needed.",
        "fixed_physics": ["flow_rate", "current", "perturbation_amplitude", "perturbation_frequency", "coulomb_softening"],
        "field_for_dynamics": "Original configured field is frozen, independently of the field mesh series.",
        "changes": "Absolute and relative changes compare adjacent successful levels; relative change denominator is the current level value.",
        "missing_results": "Missing deposition observables remain empty; failed and unrun cases remain in convergence.csv.",
        "caution": "Chaotic filament trajectories need ensemble/statistical comparisons; close integral balances alone do not establish convergence.",
        "field_only": field_only,
    }
    (out / "study_notes.json").write_text(json.dumps(instructions, indent=2) + "\n")
    executable = None
    frozen_field = None
    prerequisite_failure = None
    previous = {}
    for row, (kind, level, case) in zip(rows, designs):
        if field_only and kind != "field":
            continue
        if kind != "field" and prerequisite_failure:
            row["reason"] = "blocked_by_shared_prerequisite: " + prerequisite_failure
            _save_table(out, rows)
            continue
        directory = Path(row["directory"])
        try:
            if kind == "field":
                solve_field(case, directory)
                row.update(_field_metrics(directory, case))
                if level == "standard" and config["field_mode"] == "needle":
                    # Dynamics must retain the user's original field grid even
                    # when it differs from the standard field-study mesh.
                    if (config["field_nr"], config["field_nz"]) == (96, 160):
                        frozen_field = directory / "potential_rz.csv"
            else:
                from .cli import build, run_case
                try:
                    if executable is None:
                        executable = build(debug)
                    if config["field_mode"] == "needle":
                        if frozen_field is None:
                            frozen_field = solve_field(config, out / "frozen_dynamic_field")
                        case = {**case, "field_mode": "file", "field_file": str(frozen_field)}
                except Exception as exc:
                    # A failed shared compiler/field prerequisite is diagnosed
                    # once; do not loop identical failing setup six times.
                    prerequisite_failure = str(exc)
                    raise
                result = run_case(case, directory, executable=executable, report=False)
                row.update(_dynamic_metrics(result, directory))
            row.update(status="COMPLETED", reason="")
            _changes(row, previous.get(kind), FIELD_METRICS if kind == "field" else DYNAMIC_METRICS)
        except Exception as exc:
            row.update(status="FAILED", reason=str(exc))
        previous[kind] = row
        _save_table(out, rows)
        print(f"{kind}/{level}: {row['status']}", flush=True)
    failed = sum(row["status"] == "FAILED" for row in rows)
    completed = sum(row["status"] == "COMPLETED" for row in rows)
    print(f"{completed} completed; {failed} failed; {len(rows)-completed-failed} unrun. {out/'convergence.csv'}")
    return 1 if failed else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/demo.json")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--field-only", action="store_true", help="Run only electrostatic grid study; record dynamics as unrun")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--set", action="append", default=[], metavar="KEY=JSON")
    args = parser.parse_args(argv)
    try:
        return study(load_config(args.config, args.set), args.out, field_only=args.field_only, debug=args.debug)
    except (ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("Interrupted; convergence.csv retains completed, failed and unrun cases.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
