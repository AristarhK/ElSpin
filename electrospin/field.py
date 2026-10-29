"""Axisymmetric finite-volume Laplace field of a filled syringe needle.

The liquid meniscus is treated as equipotential.  The hollow steel tube is
filled with that liquid, so its interior does not belong to the gas domain.
This is an electrostatic boundary problem, not a Taylor-cone or charge-transport
calculation.  The curved boundary is staircased and must be mesh refined.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

EPS0 = 8.8541878128e-12


def _graded_faces(length: float, count: int, scale: float) -> np.ndarray:
    return scale * np.expm1(np.linspace(0., np.log1p(length / scale), count + 1))


def _anchor(faces: np.ndarray, value: float) -> None:
    """Align a nearby face with a known geometric interface, retaining count."""
    i = int(np.argmin(abs(faces - value)))
    if 0 < i < len(faces) - 1:
        faces[i] = value
        faces.sort()


def _solve(rf: np.ndarray, zf: np.ndarray, conductor: np.ndarray,
           voltage: float) -> dict:
    """Solve with zero side/collector potentials and insulating upstream cap.

    Conductances are face area / distance; multiplying all by epsilon0 would
    leave the solution unchanged. Symmetric elimination keeps a true SPD
    matrix. Fluxes use the identical face conductances as the matrix.
    """
    r = .5 * (rf[1:] + rf[:-1])
    z = .5 * (zf[1:] + zf[:-1])
    dr, dz = np.diff(rf), np.diff(zf)
    nr, nz = len(r), len(z)
    if conductor.shape != (nz, nr):
        raise ValueError("Conductor mask does not match grid")
    free = ~conductor
    index = np.full((nz, nr), -1, dtype=int)
    index[free] = np.arange(np.count_nonzero(free))
    n = int(np.count_nonzero(free))
    if not n or not np.any(conductor):
        raise ValueError("Need both gas cells and resolved electrode cells")
    rows, cols, vals = [], [], []
    rhs = np.zeros(n)
    faces = []
    for j, i in zip(*np.nonzero(free)):
        k = index[j, i]
        diag = 0.
        # Neighbor direction, face area, gas half-cell boundary distance.
        directions = [
            (j, i - 1, 2*np.pi*rf[i]*dz[j], .5*dr[i], -1., 0.),
            (j, i + 1, 2*np.pi*rf[i+1]*dz[j], .5*dr[i], 1., 0.),
            (j - 1, i, np.pi*(rf[i+1]**2-rf[i]**2), .5*dz[j], 0., -1.),
            (j + 1, i, np.pi*(rf[i+1]**2-rf[i]**2), .5*dz[j], 0., 1.),
        ]
        for jj, ii, area, halfdist, nr_gas, nz_gas in directions:
            if ii < 0 or jj < 0:  # symmetry axis or insulating upstream end
                continue
            external = ii >= nr or jj >= nz
            solid = not external and conductor[jj, ii]
            if external or solid:
                c = area / halfdist
                boundary_potential = voltage if solid else 0.
                diag += c
                rhs[k] += c * boundary_potential
                # Normal below is outward from electrode into the gas.
                rr = rf[i+1] if ii >= nr else (rf[i] if ii < i else r[i])
                zz = zf[j+1] if jj >= nz else (zf[j] if jj < j else z[j])
                if solid:
                    if ii != i:
                        rr, zz = (rf[max(i, ii)], z[j])
                    else:
                        rr, zz = (r[i], zf[max(j, jj)])
                electrode = "needle" if solid else ("chamber" if ii >= nr else "collector")
                faces.append((electrode, k, rr, zz, -nr_gas, -nz_gas,
                              area, halfdist, boundary_potential))
            else:
                dist = abs(r[ii]-r[i]) if ii != i else abs(z[jj]-z[j])
                c = area / dist
                diag += c
                rows.append(k); cols.append(index[jj, ii]); vals.append(-c)
        rows.append(k); cols.append(k); vals.append(diag)
    a = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    phi_free = spsolve(a, rhs)
    residual = a @ phi_free - rhs
    denominator = max(float(np.linalg.norm(rhs)), np.finfo(float).tiny)
    relative_residual = float(np.linalg.norm(residual) / denominator)
    if not np.all(np.isfinite(phi_free)) or relative_residual > 2e-9:
        raise RuntimeError(f"Field solve failed: relative residual {relative_residual:g}")
    phi = np.full((nz, nr), voltage, dtype=float)
    phi[free] = phi_free
    ez, er = np.gradient(-phi, z, r, edge_order=2)
    er[conductor] = 0.; ez[conductor] = 0.
    records = []
    charges = {"needle": 0., "collector": 0., "chamber": 0.}
    for electrode, k, rr, zz, nrr, nzz, area, distance, v in faces:
        sigma = EPS0 * (v-phi_free[k]) / distance
        charge = sigma*area
        charges[electrode] += charge
        records.append((electrode, rr, zz, nrr, nzz, area, sigma, charge))
    balance = abs(sum(charges.values())) / max(sum(abs(q) for q in charges.values()), 1e-30)
    if balance > 2e-8:
        raise RuntimeError(f"Field boundary charge balance failed: {balance:g}")
    volumes = np.diff(np.pi*rf**2)[None, :] * dz[:, None]
    return dict(r=r, z=z, phi=phi, er=er, ez=ez, conductor=conductor,
                volumes=volumes, faces=records, charges=charges,
                residual=relative_residual, charge_balance=balance)


def solve_field(config: dict, out: Path) -> Path:
    """Write needle field/charge CSVs and diagnostics; return potential CSV.

    Flat SI config keys: voltage, gap, emitter_radius (outer tube radius),
    needle_inner_radius, needle_length, field_radius, field_nr, field_nz.
    `needle_inner_radius` defaults to 0.20 mm and `emitter_radius` to 0.35 mm
    (illustrative 22G dimensions; verify the actual manufacturer's dimensions).
    """
    out = Path(out)
    voltage = float(config.get("voltage", 3500.))
    gap = float(config.get("gap", .1))
    outer = float(config.get("emitter_radius", .00035))
    inner = float(config.get("needle_inner_radius", .00020))
    length = float(config.get("needle_length", .01))
    radius = float(config.get("field_radius", .15))
    nr = int(config.get("field_nr", 96))
    nz = int(config.get("field_nz", 160))
    values = [voltage, gap, outer, inner, length, radius]
    if not all(np.isfinite(values)) or not (gap > 0 and 0 < inner < outer < radius and length > inner):
        raise ValueError("Field geometry needs 0 < inner < outer < chamber radius, gap > 0 and needle length > inner")
    if not (24 <= nr <= 600 and 48 <= nz <= 1000):
        raise ValueError("Field grid requires 24<=field_nr<=600, 48<=field_nz<=1000")
    rf = _graded_faces(radius, nr, inner / 6.)
    _anchor(rf, inner); _anchor(rf, outer)
    na = max(16, nz // 3)
    zn = -_graded_faces(length, na, inner / 6.)[::-1]
    _anchor(zn, -inner)
    zp = _graded_faces(gap, nz-na, inner / 6.)
    zf = np.concatenate([zn[:-1], zp])
    r = .5*(rf[:-1]+rf[1:]); z = .5*(zf[:-1]+zf[1:])
    rr, zz = np.meshgrid(r, z)
    # Steel annulus plus conducting liquid lumen and hemispherical meniscus.
    lower = (zz <= -inner) & (rr < outer)
    meniscus = (zz > -inner) & (zz < 0.) & (rr**2 + (zz+inner)**2 < inner**2)
    conductor = lower | meniscus
    if np.count_nonzero(meniscus) < 6:
        raise ValueError("Meniscus is unresolved; increase field_nr/field_nz")
    result = _solve(rf, zf, conductor, voltage)
    out.mkdir(parents=True, exist_ok=True)
    target = out / "potential_rz.csv"
    with target.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["r_m", "z_m", "Phi_V", "E_r_V_m", "E_z_V_m", "rho_C_m3", "cell_volume_m3"])
        for j in range(nz):
            for i in range(nr):
                writer.writerow([r[i], z[j], result["phi"][j, i], result["er"][j, i],
                                 result["ez"][j, i], 0., result["volumes"][j, i]])
    surface_charge = {"steel": 0., "liquid": 0.}
    with (out / "electrode_charge.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["electrode", "r_m", "z_m", "normal_r", "normal_z", "area_m2", "sigma_C_m2", "charge_C", "material"])
        for row in result["faces"]:
            material = "ground"
            if row[0] == "needle":
                material = "steel" if row[1] >= inner*(1-1e-10) else "liquid"
                surface_charge[material] += row[-1]
            writer.writerow([*row, material])
    conductivity = float(config.get("conductivity", 2e-6))
    permittivity = float(config.get("relative_permittivity", 40.))
    tau = EPS0*permittivity/conductivity if conductivity > 0 else None
    resolution = {
        "radial_cells_across_bore_radius": int(np.count_nonzero(r < inner)),
        "radial_cells_across_outer_radius": int(np.count_nonzero(r < outer)),
        "axial_cells_across_meniscus_height": int(np.count_nonzero((z > -inner) & (z < 0))),
        "gas_liquid_surface_facets": int(sum(row[0] == "needle" and row[1] < inner*(1-1e-10)
                                             for row in result["faces"])),
    }
    warnings = []
    if min(resolution["radial_cells_across_bore_radius"],
           resolution["axial_cells_across_meniscus_height"]) < 8:
        warnings.append("Fewer than eight cells across meniscus radius/height: refine field grid before interpreting local charge.")
    if abs(result["charges"]["collector"]) < .01*abs(result["charges"]["needle"]):
        warnings.append("Less than 1% of electrode flux reaches collector; grounded side chamber dominates this geometry.")
    summary = {
        "model": "axisymmetric Laplace; hollow steel tube filled by equipotential solution; hemispherical liquid meniscus",
        "geometry": {"outer_radius_m": outer, "inner_radius_m": inner, "needle_length_m": length,
                     "tip_z_m": 0., "metal_lip_z_m": -inner, "collector_z_m": gap,
                     "chamber_radius_m": radius},
        "boundary_conditions": {"needle_and_meniscus_V": voltage, "collector_V": 0., "chamber_V": 0.,
                                "axis": "symmetry", "upstream_outside_needle": "insulating"},
        "nr": nr, "nz": nz, "gas_cells": int(np.count_nonzero(~conductor)),
        "relative_linear_residual": result["residual"],
        "relative_boundary_charge_imbalance": result["charge_balance"],
        "electrode_charge_C": result["charges"], "needle_surface_charge_C": surface_charge,
        "min_dr_m": float(min(np.diff(rf))), "min_dz_m": float(min(np.diff(zf))),
        "geometry_resolution": resolution, "warnings": warnings,
        "electric_relaxation_time_s": tau,
        "limitations": ["Staircased metal lip and meniscus: local peak charge is mesh dependent; refine the grid.",
                        "No bevel, Taylor cone, meniscus deformation, jet or gas space charge in this field.",
                        "Equipotential liquid requires electric relaxation faster than relevant flow times.",
                        "E in CSV is cell-center reconstructed gradient; electrode charge uses conservative face fluxes."]}
    (out / "field_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return target
