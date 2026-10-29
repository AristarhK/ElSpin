"""Run the same physical problem on successive meshes. Standard library only."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import subprocess
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", default="./electrospinning")
    parser.add_argument("--config", default="config_spinning.nml")
    parser.add_argument("--out", default="mesh_study")
    parser.add_argument("--presets", nargs="+", choices=["quick", "standard", "fine", "heavy"],
                        default=["quick", "standard", "fine", "heavy"])
    args = parser.parse_args()
    exe, config = Path(args.exe).resolve(), Path(args.config).resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    results = []
    last = None
    for preset in args.presets:
        work = out / preset
        work.mkdir(exist_ok=True)
        t0 = time.perf_counter()
        proc = subprocess.run([str(exe), "--config", str(config), "--preset", preset],
                              cwd=work, text=True, capture_output=True)
        (work / "console.log").write_text(proc.stdout + proc.stderr, encoding="utf-8")
        if proc.returncode:
            raise RuntimeError(f"{preset} failed: inspect {work / 'console.log'}")
        with (work / "diagnostics.csv").open() as stream:
            d = {row["metric"]: float(row["value"]) for row in csv.DictReader(stream)}
        diameter, velocity = d["jet_diameter_end_m"], d["velocity_end_m_s"]
        row = dict(preset=preset, diameter_m=diameter, velocity_m_s=velocity,
                   temperature_C=d["temperature_end_C"], solvent_fraction=d["solvent_remaining_fraction"],
                   entry_E_V_m=d["entry_E_abs_V_m"],
                   diameter_change_percent="" if last is None else 100 * abs(diameter / last[0] - 1),
                   velocity_change_percent="" if last is None else 100 * abs(velocity / last[1] - 1),
                   min_radius_over_first_cell=d["minimum_radius_over_first_cell"],
                   poisson_residual=d["poisson_relative_residual"],
                   coupling_error=d["coupling_relative_error"], seconds=time.perf_counter()-t0)
        results.append(row)
        last = diameter, velocity
        print(f"{preset:8s} d={diameter * 1e6:.8f} um; u={velocity:.8f} m/s; "
              f"T={d['temperature_end_C']:.6f} C; time={row['seconds']:.3f} s", flush=True)
    with (out / "mesh_convergence.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=results[0])
        writer.writeheader()
        writer.writerows(results)
    if len(results)>1:
        delta = results[-1]["diameter_change_percent"]
        print(f"Final diameter change: {delta:.6g}%.")
        print("This measures grid sensitivity of this reduced model; it is not experimental accuracy.")


if __name__ == "__main__":
    main()
