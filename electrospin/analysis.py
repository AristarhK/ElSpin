"""Conservative postprocessing of computed centerlines; SI internally.

Collector maps describe bead-resolution planar chords. They do not resolve fiber
contacts, pore sizes, stacking, or sub-bead coiling. Historical parcel charges
are charges delivered to ground, not residual charge stored on deposited fibers.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np

EPSILON_0 = 8.8541878128e-12


def read_table(path: Path) -> np.ndarray:
    """Read named numerical CSV, including a valid header-only result."""
    path = Path(path)
    if not path.exists():
        return np.empty(0, dtype=[])
    with path.open(encoding="utf-8") as stream:
        reader = csv.reader(stream)
        try:
            names = [name.strip() for name in next(reader)]
        except StopIteration:
            return np.empty(0, dtype=[])
        values = []
        for line_number, row in enumerate(reader, 2):
            if not row or not any(value.strip() for value in row):
                continue
            if len(row) != len(names):
                raise ValueError(f"{path}:{line_number}: expected {len(names)} CSV fields, found {len(row)}")
            try:
                values.append([float(value.replace("D", "E")) for value in row])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid numeric CSV field") from exc
    result = np.empty(len(values), dtype=[(name, float) for name in names])
    for j, name in enumerate(names):
        result[name] = [row[j] for row in values]
    if any(not np.isfinite(result[name]).all() for name in names):
        raise ValueError(f"Non-finite data in {path}")
    return result


def deposited_segments(deposition: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Consecutive deposited IDs only; length and material mass belong to i→i+1.

    The core stores segment i between parcel i and i+1. An isolated hit does not
    establish a fiber chord. Sorting by impact time would create false fibers.
    """
    if len(deposition) < 2:
        return np.empty((0, 2, 2)), np.empty(0), np.empty(0)
    data = np.sort(deposition, order="bead_id")
    valid = np.diff(data["bead_id"]) == 1
    xy = np.column_stack([data["collector_x_m"], data["collector_y_m"]])
    segments = np.stack([xy[:-1][valid], xy[1:][valid]], axis=1)
    masses = data["polymer_mass_kg"][:-1][valid]
    lengths = np.linalg.norm(segments[:, 1] - segments[:, 0], axis=1)
    return segments, masses, lengths


def coverage_map(deposition: np.ndarray, bins: int = 64) -> dict:
    """Rasterize chord length and polymer mass with conservative midpoint weights.

    Per-segment quadrature has ≤ half-cell spacing. Every sample contributes its
    share of segment mass and length. Cell density integrates back to those sums.
    Grid extent follows the predicted material-frame impact region, not a claimed
    collector boundary. Coincident points retain mass but carry zero line length.
    """
    if bins < 4:
        raise ValueError("coverage bins must be >= 4")
    segments, masses, lengths = deposited_segments(deposition)
    if len(deposition):
        points = np.column_stack([deposition["collector_x_m"], deposition["collector_y_m"]])
        center = (points.max(axis=0) + points.min(axis=0)) / 2
        span = max(float(np.ptp(points, axis=0).max()), 1e-5) * 1.1
    else:
        center, span = np.zeros(2), .002
    xedges = np.linspace(center[0] - span / 2, center[0] + span / 2, bins + 1)
    yedges = np.linspace(center[1] - span / 2, center[1] + span / 2, bins + 1)
    mass_grid = np.zeros((bins, bins))
    length_grid = np.zeros_like(mass_grid)
    dx = span / bins
    for segment, mass, length in zip(segments, masses, lengths):
        count = max(1, int(np.ceil(length / (0.5 * dx))))
        fractions = (np.arange(count) + .5) / count
        points = segment[0] + fractions[:, None] * (segment[1] - segment[0])
        ix = np.clip(np.searchsorted(xedges, points[:, 0], side="right") - 1, 0, bins - 1)
        iy = np.clip(np.searchsorted(yedges, points[:, 1], side="right") - 1, 0, bins - 1)
        np.add.at(mass_grid, (iy, ix), mass / count)
        np.add.at(length_grid, (iy, ix), length / count)
    return {"xedges": xedges, "yedges": yedges, "cell_area_m2": dx * dx,
            "mass_kg": mass_grid, "length_m": length_grid,
            "areal_polymer_mass_kg_m2": mass_grid / (dx * dx),
            "line_density_m_m2": length_grid / (dx * dx)}


def field_data(run_dir: Path) -> np.ndarray:
    for name in ("field/potential_rz.csv", "potential_rz.csv"):
        data = read_table(Path(run_dir) / name)
        if len(data):
            return data
    return np.empty(0, dtype=[])


def electrode_charge(run_dir: Path, electrode: str = "needle") -> np.ndarray:
    """Select named electrode faces from the optional geometry-aware field solver."""
    path = Path(run_dir) / "field/electrode_charge.csv"
    names = ("r_m", "z_m", "normal_r", "normal_z", "area_m2", "sigma_C_m2", "charge_C")
    if not path.exists():
        return np.empty(0, dtype=[(name, float) for name in names])
    with path.open(encoding="utf-8") as stream:
        rows = [row for row in csv.DictReader(stream) if row["electrode"] == electrode]
    result = np.empty(len(rows), dtype=[(name, float) for name in names])
    for name in names:
        result[name] = [float(row[name]) for row in rows]
    return result


def field_grid(data: np.ndarray, name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Recover rectangular original axisymmetric grid without row-order assumptions."""
    r, z = np.unique(data["r_m"]), np.unique(data["z_m"])
    if len(r) * len(z) != len(data):
        raise ValueError("Field table is not a complete rectangular r-z grid")
    values = np.full((len(z), len(r)), np.nan)
    values[np.searchsorted(z, data["z_m"]), np.searchsorted(r, data["r_m"])] = data[name]
    if not np.isfinite(values).all():
        raise ValueError("Missing or duplicated field grid points")
    return r, z, values


def _write_coverage(run_dir: Path, coverage: dict) -> None:
    xs = (coverage["xedges"][:-1] + coverage["xedges"][1:]) / 2
    ys = (coverage["yedges"][:-1] + coverage["yedges"][1:]) / 2
    with (run_dir / "coverage.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["collector_x_m", "collector_y_m", "cell_area_m2",
                         "polymer_areal_mass_kg_m2", "line_density_m_m2"])
        for iy, y in enumerate(ys):
            for ix, x in enumerate(xs):
                writer.writerow([x, y, coverage["cell_area_m2"],
                                 coverage["areal_polymer_mass_kg_m2"][iy, ix],
                                 coverage["line_density_m_m2"][iy, ix]])


def _write_vtp(run_dir: Path, trajectory: np.ndarray) -> None:
    """Export final bead chain and point fields in standard VTK PolyData XML."""
    if not len(trajectory):
        return
    final = trajectory[trajectory["time_s"] == trajectory["time_s"].max()]
    final = np.sort(final, order="bead_id")
    edges = np.flatnonzero(np.diff(final["bead_id"]) == 1)
    vtk = ET.Element("VTKFile", type="PolyData", version="0.1", byte_order="LittleEndian")
    piece = ET.SubElement(ET.SubElement(vtk, "PolyData"), "Piece", NumberOfPoints=str(len(final)),
                          NumberOfVerts="0", NumberOfLines=str(len(edges)),
                          NumberOfStrips="0", NumberOfPolys="0")

    def array(parent, name, values, dtype="Float64", components=1):
        element = ET.SubElement(parent, "DataArray", type=dtype, Name=name,
                                NumberOfComponents=str(components), format="ascii")
        element.text = " ".join(str(v) for v in np.asarray(values).ravel())

    coordinates = np.column_stack([final[name] for name in ("x_m", "y_m", "z_m")])
    array(ET.SubElement(piece, "Points"), "Points", coordinates, components=3)
    lines = ET.SubElement(piece, "Lines")
    array(lines, "connectivity", np.column_stack([edges, edges + 1]), "Int64")
    array(lines, "offsets", 2 * (np.arange(len(edges)) + 1), "Int64")
    point_data = ET.SubElement(piece, "PointData", Scalars="radius_m")
    for name in ("bead_id", "radius_m", "charge_C", "temperature_K", "stress_Pa", "deposited"):
        array(point_data, name, final[name])
    ET.ElementTree(vtk).write(run_dir / "trajectory.vtp", encoding="utf-8", xml_declaration=True)


def _write_final_jet(run_dir: Path, trajectory: np.ndarray) -> None:
    """Export segment-average diagnostics; no extra electrostatic closure implied."""
    if not len(trajectory) or "segment_length_m" not in trajectory.dtype.names:
        return
    final = np.sort(trajectory[trajectory["time_s"] == trajectory["time_s"].max()], order="bead_id")
    names = ["bead_id", "x_m", "y_m", "z_m", "segment_length_m", "wet_diameter_m",
             "speed_m_s", "polymer_mass_fraction", "linear_charge_C_m",
             "surface_charge_segment_average_C_m2", "temperature_K", "stress_Pa", "deposited"]
    with (run_dir / "final_jet.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(names)
        for row in final:
            length = row["segment_length_m"]
            radius = row["radius_m"]
            linear = row["charge_C"] / length if length > 0 else None
            surface = linear / (2 * np.pi * radius) if linear is not None and radius > 0 else None
            speed = np.sqrt(sum(row[name] ** 2 for name in ("vx_m_s", "vy_m_s", "vz_m_s")))
            fraction = row["polymer_mass_kg"] / row["mass_kg"] if row["mass_kg"] > 0 else None
            writer.writerow([row["bead_id"], row["x_m"], row["y_m"], row["z_m"], length,
                             radius * 2, speed, fraction, linear, surface, row["temperature_K"],
                             row["stress_Pa"], row["deposited"]])


def summarize_run(run_dir: Path) -> dict:
    """Write reproducible analysis.json, coverage.csv, emitter_charge.csv and VTP."""
    run_dir = Path(run_dir)
    trajectory = read_table(run_dir / "trajectory.csv")
    deposition = read_table(run_dir / "deposition.csv")
    diagnostics = read_table(run_dir / "diagnostics.csv")
    config_path = run_dir / "config.json"
    config = json.loads(config_path.read_text(encoding="utf-8")) if config_path.exists() else {}
    segments, masses, lengths = deposited_segments(deposition)
    coverage = coverage_map(deposition)
    _write_coverage(run_dir, coverage)
    _write_vtp(run_dir, trajectory)
    _write_final_jet(run_dir, trajectory)
    nonzero = lengths > 0
    orientation = None
    mean_angle = None
    if np.any(nonzero):
        delta = segments[nonzero, 1] - segments[nonzero, 0]
        angles = np.arctan2(delta[:, 1], delta[:, 0])
        order = np.sum(lengths[nonzero] * np.exp(2j * angles)) / np.sum(lengths[nonzero])
        orientation = float(abs(order))
        mean_angle = float(np.angle(order) * 90 / np.pi) if orientation > 1e-8 else None
    result = {
        "deposited_parcels": len(deposition),
        "connected_deposited_chords": len(segments),
        "mapped_polymer_mass_kg": float(masses.sum()),
        "delivered_polymer_mass_kg": float(deposition["polymer_mass_kg"].sum()) if len(deposition) else 0.,
        "mapped_centerline_length_m": float(lengths.sum()),
        "orientation_nematic_order_0_to_1": orientation,
        "principal_orientation_deg": mean_angle,
        "mean_wet_diameter_at_hits_m": float(2 * np.mean(deposition["radius_m"])) if len(deposition) else None,
        "mean_dry_diameter_same_length_at_hits_m": float(2 * np.mean(deposition["dry_radius_m"])) if len(deposition) else None,
        "delivered_charge_C": float(deposition["charge_C"].sum()) if len(deposition) else 0.,
        "coverage_mass_integral_error_kg": float(coverage["mass_kg"].sum() - masses.sum()),
        "notes": [
            "Collector map uses consecutive deposited bead centerline chords in material coordinates.",
            "Unconnected landing parcels have no resolved chord and are excluded from the coverage map.",
            "Chord polymer mass belongs to the lower bead ID; no fiber contact, stacking or pore-size model.",
            "Dry diameter assumes solvent removal at unchanged local segment length.",
            "Delivered charge is transferred to ideal ground; it is not a retained surface-charge map.",
            "Hit diameter averages are parcel-count weighted, not experimentally calibrated distributions.",
            "final_jet.csv charge densities assume parcel charge uniformly distributed over its segment length and cylindrical surface; deposited values are historical delivered charges.",
            "3D chain diameter colors show current segment-equivalent diameters; impact histograms use diameters recorded at crossing, and static 3D transverse axes are expanded.",
        ],
        "applicability_warnings": [
            "3D convective parcel charge is prescribed; a full leaky-dielectric charge-transport law and non-axisymmetric electrical feedback are not solved.",
        ],
    }
    if len(deposition):
        result["first_impact_time_s"] = float(deposition["time_s"].min())
        radial = np.hypot(deposition["collector_x_m"], deposition["collector_y_m"])
        result["rms_impact_radius_material_m"] = float(np.sqrt(np.mean(radial ** 2)))
        result["max_impact_radius_material_m"] = float(radial.max())
    if len(trajectory):
        result["last_snapshot_time_s"] = float(trajectory["time_s"].max())
        final = trajectory[trajectory["time_s"] == trajectory["time_s"].max()]
        live = final[final["deposited"] < .5]
        result["final_active_beads"] = len(live)
        result["final_max_radial_excursion_m"] = float(np.hypot(live["x_m"], live["y_m"]).max()) if len(live) else 0.
        if len(live) and all(name in live.dtype.names for name in ("vx_m_s", "vy_m_s", "vz_m_s")):
            speed = np.sqrt(live["vx_m_s"] ** 2 + live["vy_m_s"] ** 2 + live["vz_m_s"] ** 2)
            result["final_active_mean_speed_m_s"] = float(speed.mean())
            result["final_active_max_speed_m_s"] = float(speed.max())
            result["final_active_mean_polymer_mass_fraction"] = float(np.mean(live["polymer_mass_kg"] / live["mass_kg"]))
    if float(config.get("conductivity", 0)) > 0:
        result["input_electrical_relaxation_screening_time_s"] = EPSILON_0 * float(config.get("relative_permittivity", 1)) / float(config["conductivity"])
        result["notes"].append("Electrical relaxation time epsilon0*epsilon_r/K is an input-property screening diagnostic; dielectric polarization is not resolved.")
    if config.get('charge_transport',False):
        result['applicability_warnings']=[w for w in result['applicability_warnings'] if not w.startswith('3D convective parcel charge')]
        result['applicability_warnings'].append('Reduced longitudinal Ohmic transport uses positive-definite regularized bead capacitances, prescribed inlet total current and grounded deposited nodes. Resolve injection interval and softening; this is not full leaky-dielectric transport.')
    if all(name in config for name in ("flow_rate", "jet_radius", "density_polymer", "density_solvent", "polymer_fraction")):
        fraction = float(config["polymer_fraction"])
        density = 1 / (fraction / float(config["density_polymer"]) + (1 - fraction) / float(config["density_solvent"]))
        result["inlet_mean_speed_m_s"] = float(config["flow_rate"]) / (np.pi * float(config["jet_radius"]) ** 2)
        result["inlet_mass_flow_kg_s"] = density * float(config["flow_rate"])
        result["inlet_polymer_mass_flow_kg_s"] = result["inlet_mass_flow_kg_s"] * fraction
    duration = float(result.get("last_snapshot_time_s", 0))
    if duration > 0:
        result["time_averaged_delivered_polymer_rate_kg_s"] = result["delivered_polymer_mass_kg"] / duration
        result["time_averaged_delivered_convective_current_A"] = result["delivered_charge_C"] / duration
        result["notes"].append("Delivered rates average the finite startup simulation window and are not steady-state collection efficiencies.")
    if len(diagnostics):
        result["last_diagnostics"] = {name: float(diagnostics[-1][name]) for name in diagnostics.dtype.names}
        for name in ("mass_relative_error", "polymer_relative_error", "charge_relative_error"):
            if name in diagnostics.dtype.names:
                result[f"maximum_{name}"] = float(np.max(np.abs(diagnostics[name])))
    field = field_data(run_dir)
    if len(field):
        needle = electrode_charge(run_dir)
        if len(needle):
            result["needle_surface_charge_C"] = float(needle["charge_C"].sum())
            result["notes"].append("Needle surface charge uses finite-volume electrode-face flux on a staircase approximation; geometric mesh convergence is required.")
        else:
            first = np.sort(field[field["z_m"] == field["z_m"].min()], order="r_m")
            first = first[first["r_m"] <= float(config.get("emitter_radius", 2.5e-4))]
            with (run_dir / "emitter_charge.csv").open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(["r_m", "sample_z_m", "E_z_first_cell_V_m", "sigma_disk_first_cell_estimate_C_m2"])
                for row in first:
                    writer.writerow([row["r_m"], row["z_m"], row["E_z_V_m"], EPSILON_0 * row["E_z_V_m"]])
            result["notes"].append("Emitter charge is epsilon0*Ez in the first gas cell under an emitting disk, not resolved needle surface charge; refine the field mesh.")
        result["field_max_magnitude_V_m"] = float(np.hypot(field["E_r_V_m"], field["E_z_V_m"]).max())
        peak = field[int(np.argmax(np.hypot(field["E_r_V_m"], field["E_z_V_m"])))]
        result["field_peak_location_m"] = {"r": float(peak["r_m"]), "z": float(peak["z_m"])}
        result["applicability_warnings"].extend([
            "Peak field and charge near a sharp electrode lip depend on staircase geometry and mesh resolution; a single-grid maximum is not a converged tip-field prediction.",
            "Gas discharge is not coupled to the 3D jet/field calculation. High local fields require a separate corona assessment; these plots do not establish discharge onset.",
        ])
    (run_dir / "analysis.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return result
