"""Revolve the *computed* axisymmetric electrode facets, preserving their sigma.

The exported surface is a 3D rendering of a 2D r-z electrostatic solution. It is
not a separate 3D Poisson solve, CAD refinement, or fitted Taylor-cone shape.
"""
from __future__ import annotations

import csv
import html
import json
import os
from pathlib import Path
import tempfile

import numpy as np

from .analysis import field_data, field_grid


def _faces(run_dir: Path) -> list[dict]:
    source = run_dir / "field/electrode_charge.csv"
    if not source.exists():
        return []
    with source.open(encoding="utf-8") as stream:
        rows = [row for row in csv.DictReader(stream) if row["electrode"] == "needle"]
    result = []
    for row in rows:
        result.append({name: row[name] if name in ("electrode", "material") else float(row[name]) for name in row})
        result[-1].setdefault("material", "electrode")
    return result


def facet_mesh(faces: list[dict], sectors: int = 40) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return points, triangles and per-triangle sigma of revolved FV faces.

    r is the arithmetic midpoint on an annular face, so area=2*pi*r*dr
    reconstructs its exact radial width. Cylindrical faces use area=2*pi*r*dz.
    Polygonal angular subdivision approximates the ring; it adds no new field
    degrees of freedom. Input charge_C/area stays constant on each FV face.
    """
    if sectors < 8:
        raise ValueError("At least eight angular sectors are required")
    points, triangles, values = [], [], []
    angles = np.linspace(0, 2 * np.pi, sectors + 1)
    for face in faces:
        r, z, area = face["r_m"], face["z_m"], face["area_m2"]
        if r <= 0 or area <= 0:
            raise ValueError("Electrode facet has non-positive radius or area")
        width = area / (2 * np.pi * r)
        if abs(face["normal_r"]) > .5:
            radius = np.array([r, r])
            axial = np.array([z - width / 2, z + width / 2])
        elif abs(face["normal_z"]) > .5:
            radius = np.maximum([r - width / 2, r + width / 2], 0.)
            axial = np.array([z, z])
        else:
            raise ValueError("Only radial/axial finite-volume facets are supported")
        offset = len(points)
        for rr, zz in zip(radius, axial):
            points.extend(np.column_stack([rr * np.cos(angles), rr * np.sin(angles), np.full(len(angles), zz)]))
        stride = sectors + 1
        for j in range(sectors):
            triangles.extend([[offset + j, offset + j + 1, offset + stride + j],
                              [offset + j + 1, offset + stride + j + 1, offset + stride + j]])
            values.extend([face["sigma_C_m2"]] * 2)
    return np.asarray(points).reshape(-1, 3), np.asarray(triangles, dtype=int).reshape(-1, 3), np.asarray(values)


def render_needle(run_dir: Path) -> list[Path]:
    """Create independent self-contained needle3d.html and needle3d.png."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    run_dir = Path(run_dir)
    faces = _faces(run_dir)
    if not faces:
        return []
    summary_path = run_dir / "field/field_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.exists() else {}
    geometry = summary.get("geometry", {})
    outer = float(geometry.get("outer_radius_m", max(face["r_m"] for face in faces)))
    tip = float(geometry.get("tip_z_m", 0))
    near = [face for face in faces if face["z_m"] >= tip - 3 * outer]
    if not near:
        near = faces
    groups = sorted({face["material"] for face in faces})
    labels = {"steel": "Steel tube / lip", "liquid": "Conducting liquid meniscus", "electrode": "Electrode"}
    all_sigma = np.array([face["sigma_C_m2"] for face in faces]) * 1e6
    cmin, cmax = float(all_sigma.min()), float(all_sigma.max())
    if cmax == cmin:
        cmax = cmin + max(abs(cmin) * .1, 1e-12)
    figure = make_subplots(rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "scene"}]],
                           subplot_titles=("Computed electrode facets · full needle", "Near-tip steel lip and liquid meniscus"))
    for col, subset in enumerate((faces, near), 1):
        for material in groups:
            selected = [face for face in subset if face["material"] == material]
            if not selected:
                continue
            points, triangles, sigma = facet_mesh(selected)
            figure.add_trace(go.Mesh3d(
                x=points[:, 0] * 1e3, y=points[:, 1] * 1e3, z=points[:, 2] * 1e3,
                i=triangles[:, 0], j=triangles[:, 1], k=triangles[:, 2],
                intensity=sigma * 1e6, intensitymode="cell", cmin=cmin, cmax=cmax,
                colorscale="Viridis", showscale=col == 2 and material == groups[-1],
                colorbar=dict(title="σ [µC/m²]", len=.72), flatshading=True,
                lighting=dict(ambient=.82, diffuse=.6, specular=.1, roughness=.9),
                name=labels.get(material, material), legendgroup=material, showlegend=col == 2,
                text=labels.get(material, material),
                hovertemplate="%{text}<br>x=%{x:.5g} mm<br>y=%{y:.5g} mm<br>z=%{z:.5g} mm<br>σ=%{intensity:.5g} µC/m²<extra></extra>"), row=1, col=col)
    scene = dict(xaxis_title="x [mm]", yaxis_title="y [mm]", zaxis_title="z [mm]",
                 aspectmode="data", zaxis=dict(autorange="reversed"),
                 camera=dict(eye=dict(x=1.6, y=1.5, z=.6)))
    figure.update_layout(template="plotly_white", height=720, scene=scene, scene2=scene,
                         legend=dict(orientation="h", y=-.05), margin=dict(t=65, b=80))
    notes = [
        "Поверхность получена вращением рассчитанных осесимметричных граней, а не отдельным решением 3D-Пуассона.",
        "Цвет показывает заряд на конечном объёме: σ=ε₀Eₙ, вычисленный консервативным потоком на грани проводника.",
        "Сталь и проводящая жидкость имеют общий потенциал. Форма мениска задана полусферой; конус Тейлора и деформация поверхности не рассчитываются.",
        "Ступенчатая геометрия губки и мениска сохранена; вращение не улучшает разрешение поля. Пиковый заряд и поле у острых краёв чувствительны к сетке.",
        "Газовый разряд не связан с этой моделью поля. Высокое локальное поле требует отдельной оценки короны; карта не доказывает возникновение разряда.",
    ]
    notes_html = "".join(f"<li>{html.escape(note)}</li>" for note in notes)
    document = """<!doctype html><html lang="ru"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Needle surface charge · axisymmetric reconstruction</title><style>body{font-family:system-ui,sans-serif;background:#f4f6fc;color:#172033;margin:0}main{max-width:1450px;padding:28px;margin:auto}li,p{line-height:1.6}.plot{background:white;border-radius:16px}h1{margin-bottom:8px}</style><main>
<h1>Заряд на игле и жидком мениске</h1><p>Интерактивная 3D-реконструкция осесимметричного электростатического расчёта. Масштаб осей геометрически согласован. В легенде можно скрыть сталь или жидкость.</p><div class="plot">""" + figure.to_html(full_html=False, include_plotlyjs=True, config={"responsive": True, "displaylogo": False}) + \
        "</div><ul>" + notes_html + "</ul><p>Исходные значения: field/electrode_charge.csv; геометрия и невязки: field/field_summary.json.</p></main></html>"
    html_path = run_dir / "needle3d.html"
    html_path.write_text(document, encoding="utf-8")

    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "electrospin-mpl"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    fig = plt.figure(figsize=(15, 6), layout="constrained", facecolor="#f4f6fc")
    spec = fig.add_gridspec(1, 3, width_ratios=[.9, 1.1, 1.1])
    norm = Normalize(cmin, cmax)
    for column, subset in enumerate((faces, near)):
        ax = fig.add_subplot(spec[0, column], projection="3d")
        points, triangles, sigma = facet_mesh(subset, sectors=24)
        polygons = points[triangles] * 1e3
        artist = Poly3DCollection(polygons, facecolors=plt.get_cmap("viridis")(norm(sigma * 1e6)),
                                  edgecolors="none", linewidths=0, antialiased=False)
        ax.add_collection3d(artist)
        bounds = points * 1e3
        for setter, axis in ((ax.set_xlim, 0), (ax.set_ylim, 1), (ax.set_zlim, 2)):
            lo, hi = bounds[:, axis].min(), bounds[:, axis].max()
            pad = max((hi - lo) * .02, 1e-6)
            setter(lo - pad, hi + pad)
        ax.set_box_aspect(np.maximum(np.ptp(bounds, axis=0), 1e-8))
        ax.invert_zaxis()
        ax.view_init(elev=12, azim=-55)
        ax.set(xlabel="x [mm]", ylabel="y [mm]", zlabel="z [mm]",
               title="Full needle" if column == 0 else "Computed near-tip facets")
        if column == 0:
            # Equal geometric scaling makes transverse ticks overlap on a long tube.
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.set_title(f"Full needle · Ø{outer * 2e3:.3g} mm")
        ax.tick_params(labelsize=7)
    ax = fig.add_subplot(spec[0, 2])
    field = field_data(run_dir)
    if len(field):
        r, z, er = field_grid(field, "E_r_V_m")
        _, _, ez = field_grid(field, "E_z_V_m")
        ir = r <= 3 * outer
        iz = (z >= tip - 3 * outer) & (z <= tip + 3 * outer)
        if ir.any() and iz.any():
            mag = np.hypot(er, ez)[np.ix_(iz, ir)]
            mirrored = np.column_stack([mag[:, ::-1], mag])
            loge = np.ma.log10(np.ma.masked_less_equal(mirrored, 0))
            xs = np.r_[-r[ir][::-1], r[ir]] * 1e3
            image = ax.pcolormesh(xs, z[iz] * 1e3, loge, shading="nearest", cmap="magma")
            fig.colorbar(image, ax=ax, label="log10 |E| [V/m]", shrink=.8)
            ax.invert_yaxis()
    ax.set(xlabel="Meridional x [mm]", ylabel="z [mm]", title="Near-tip gas electric field")
    ax.set_aspect("equal")
    scalar = ScalarMappable(norm=norm, cmap="viridis")
    fig.colorbar(scalar, ax=fig.axes[:2], label="Surface charge σ [µC/m²]", shrink=.65, pad=.03)
    fig.suptitle("Needle charge · revolved axisymmetric FV solution (not a 3D field solve)", fontsize=14)
    png_path = run_dir / "needle3d.png"
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    return [html_path, png_path]
