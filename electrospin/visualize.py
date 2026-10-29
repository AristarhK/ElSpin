"""Offline numerical reports for electrospinning runs (no synthetic fiber paths)."""
from __future__ import annotations

import html
import json
import os
from pathlib import Path
import tempfile

import numpy as np

from .analysis import (coverage_map, deposited_segments, electrode_charge, field_data,
                       field_grid, read_table, summarize_run)


def _chain(data):
    """Insert gaps where IDs cease to be consecutive."""
    data = np.sort(data, order="bead_id")
    x, y, z = [], [], []
    previous = None
    for row in data:
        if previous is not None and row["bead_id"] != previous + 1:
            x.append(None); y.append(None); z.append(None)
        x.append(float(row["x_m"]) * 1e3)
        y.append(float(row["y_m"]) * 1e3)
        z.append(float(row["z_m"]) * 1e3)
        previous = row["bead_id"]
    return x, y, z


def _line_xy(segments):
    xs, ys = [], []
    for segment in segments:
        xs.extend([segment[0, 0] * 1e3, segment[1, 0] * 1e3, None])
        ys.extend([segment[0, 1] * 1e3, segment[1, 1] * 1e3, None])
    return xs, ys


def _static_report(run_dir, trajectory, deposits, field, coverage):
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "electrospin-mpl"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titlesize": 11,
                         "axes.labelsize": 10, "figure.facecolor": "#f6f7fb"})
    fig = plt.figure(figsize=(17, 11), layout="constrained")
    grid = fig.add_gridspec(2, 3)
    ax = fig.add_subplot(grid[0, 0], projection="3d")
    if len(trajectory):
        times = np.unique(trajectory["time_s"])
        picks = np.unique(np.linspace(0, len(times) - 1, min(6, len(times))).astype(int))
        for index in picks[:-1]:
            row = trajectory[trajectory["time_s"] == times[index]]
            x, y, z = _chain(row)
            ax.plot(x, y, z, color="#9fa8da", alpha=.25, linewidth=.7)
        last = trajectory[trajectory["time_s"] == times[-1]]
        x, y, z = _chain(last)
        ax.plot(x, y, z, color="#4338ca", linewidth=1)
        scatter = ax.scatter(last["x_m"] * 1e3, last["y_m"] * 1e3,
                             last["z_m"] * 1e3, c=2e6 * last["radius_m"],
                             cmap="plasma", s=9)
        fig.colorbar(scatter, ax=ax, label="Segment equivalent diameter [µm]", shrink=.6)
    ax.set(xlabel="x [mm]", ylabel="y [mm]", zlabel="z [mm]", title="3D chain · transverse axes expanded")
    ax.invert_zaxis()

    ax = fig.add_subplot(grid[0, 1])
    segments, masses, lengths = deposited_segments(deposits)
    if len(segments):
        ax.add_collection(LineCollection(segments * 1e3, colors="#4338ca", linewidths=.9))
        ax.autoscale()
    if len(deposits):
        ax.scatter(deposits["collector_x_m"] * 1e3, deposits["collector_y_m"] * 1e3,
                   s=6, c="#f97316", label="Impact parcels")
        ax.legend(fontsize=8)
    else:
        ax.text(.5, .5, "No collector impacts in simulated time", ha="center", transform=ax.transAxes)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
    ax.set(xlabel="Material x [mm]", ylabel="Material y [mm]", title="Collector: consecutive deposited chords")
    ax.set_aspect("equal", adjustable="datalim")

    ax = fig.add_subplot(grid[0, 2])
    image = ax.pcolormesh(coverage["xedges"] * 1e3, coverage["yedges"] * 1e3,
                         coverage["areal_polymer_mass_kg_m2"] * 1e6, cmap="magma", shading="flat",
                         vmin=0, vmax=max(float(coverage["areal_polymer_mass_kg_m2"].max() * 1e6), 1e-30))
    fig.colorbar(image, ax=ax, label="Polymer areal mass [mg/m²]", shrink=.8)
    ax.set(xlabel="Material x [mm]", ylabel="Material y [mm]", title="Resolved deposition mass / cell area")
    ax.set_aspect("equal")

    ax = fig.add_subplot(grid[1, 0])
    if len(deposits):
        ax.hist(2e6 * deposits["radius_m"], bins=min(25, max(4, int(np.sqrt(len(deposits))))),
                color="#6366f1", alpha=.75, label="Wet at impact")
        ax.hist(2e6 * deposits["dry_radius_m"], bins=min(25, max(4, int(np.sqrt(len(deposits))))),
                histtype="step", color="#f97316", linewidth=1.5, label="Dry, unchanged length")
        ax.legend(fontsize=8)
    ax.set(xlabel="Diameter [µm]", ylabel="Parcel count", title="Predicted impact diameter distribution")

    ax = fig.add_subplot(grid[1, 1])
    if len(field):
        r, z, er = field_grid(field, "E_r_V_m")
        _, _, ez = field_grid(field, "E_z_V_m")
        magnitude = np.hypot(er, ez)
        image = ax.pcolormesh(r * 1e3, z * 1e3, np.ma.log10(np.ma.masked_less_equal(magnitude, 0)),
                             cmap="viridis", shading="nearest")
        fig.colorbar(image, ax=ax, label="log10 |E| [V/m]", shrink=.8)
        ax.invert_yaxis()
    else:
        ax.text(.5, .5, "Uniform-field model: no field mesh", ha="center", transform=ax.transAxes)
    ax.set(xlabel="r [mm]", ylabel="z [mm]", title="Axisymmetric external electric field")

    ax = fig.add_subplot(grid[1, 2])
    needle = electrode_charge(run_dir)
    if len(needle):
        order = np.argsort(needle["z_m"])
        ax.scatter(needle["z_m"][order] * 1e3, needle["sigma_C_m2"][order] * 1e6,
                   c=needle["r_m"][order] * 1e3, cmap="viridis", s=10)
        ax.set(xlabel="Needle face z [mm]", ylabel="Surface charge [µC/m²]",
               title="Needle face charge: staircase geometry")
    else:
        proxy = read_table(run_dir / "emitter_charge.csv")
        if len(proxy):
            ax.plot(proxy["r_m"] * 1e3, proxy["sigma_disk_first_cell_estimate_C_m2"] * 1e6,
                    color="#4338ca")
        ax.set(xlabel="Disk radius [mm]", ylabel="First-cell ε₀Ez [µC/m²]",
               title="Disk normal-charge proxy; mesh sensitive")
    fig.suptitle("Electrospinning · computed trajectories, impact chords and field", fontsize=17, fontweight="bold")
    path = run_dir / "report.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def render_report(run_dir: Path) -> list[Path]:
    """Generate self-contained report.html, report.png and numerical exports."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    run_dir = Path(run_dir)
    summary = summarize_run(run_dir)
    trajectory = read_table(run_dir / "trajectory.csv")
    deposits = read_table(run_dir / "deposition.csv")
    diagnostics = read_table(run_dir / "diagnostics.csv")
    field = field_data(run_dir)
    coverage = coverage_map(deposits)
    config = json.loads((run_dir / "config.json").read_text(encoding="utf-8")) if (run_dir / "config.json").exists() else {}
    template = "plotly_white"
    figures = []

    flight = go.Figure()
    if len(trajectory):
        times = np.unique(trajectory["time_s"])
        selection = np.unique(np.linspace(0, len(times) - 1, min(60, len(times))).astype(int))

        def trace(time):
            rows = trajectory[trajectory["time_s"] == time]
            x, y, z = _chain(rows)
            return go.Scatter3d(x=x, y=y, z=z, mode="lines", line=dict(color="#4f46e5", width=4),
                                name="Connected bead chain", hoverinfo="x+y+z")

        def beads(time):
            rows = trajectory[trajectory["time_s"] == time]
            return go.Scatter3d(x=rows["x_m"] * 1e3, y=rows["y_m"] * 1e3, z=rows["z_m"] * 1e3,
                                mode="markers", marker=dict(size=3, color=rows["radius_m"] * 2e6,
                                colorscale="Plasma", cmin=0, cmax=float(trajectory["radius_m"].max() * 2e6),
                                colorbar=dict(title="Segment d [µm]")),
                                customdata=np.column_stack([rows["bead_id"], rows["temperature_K"], rows["stress_Pa"]]),
                                hovertemplate="Bead %{customdata[0]:.0f}<br>x=%{x:.4g} mm<br>y=%{y:.4g} mm<br>z=%{z:.4g} mm<br>T=%{customdata[1]:.3f} K<br>stress=%{customdata[2]:.3g} Pa<extra></extra>",
                                name="Computed parcels")

        flight.add_trace(trace(times[-1]))
        flight.add_trace(beads(times[-1]))
        frames = [go.Frame(name=str(j), data=[trace(times[i]), beads(times[i])], traces=[0, 1])
                  for j, i in enumerate(selection)]
        flight.frames = frames
        flight.update_layout(sliders=[dict(active=len(selection) - 1, currentvalue=dict(prefix="t [s] = "),
                              steps=[dict(label=f"{times[i]:.5g}", method="animate", args=[[str(j)],
                                     dict(mode="immediate", frame=dict(duration=0, redraw=True), transition=dict(duration=0))])
                                     for j, i in enumerate(selection)])],
                             updatemenus=[dict(type="buttons", buttons=[
                                 dict(label="▶ Play", method="animate", args=[None, dict(fromcurrent=False, frame=dict(duration=100, redraw=True), transition=dict(duration=0))]),
                                 dict(label="Pause", method="animate", args=[[None], dict(mode="immediate", frame=dict(duration=0, redraw=False))])])])
    gap = float(config.get("gap", .1))
    extent = max(float(np.max(np.abs(trajectory["x_m"]))) if len(trajectory) else .001,
                 float(np.max(np.abs(trajectory["y_m"]))) if len(trajectory) else .001, .001) * 1.2e3
    flight.add_trace(go.Mesh3d(x=[-extent, extent, extent, -extent], y=[-extent, -extent, extent, extent],
                               z=[gap * 1e3] * 4, i=[0, 0], j=[1, 2], k=[2, 3],
                               opacity=.18, color="#94a3b8", name="Collector z=L", hoverinfo="name"))
    flight.update_layout(template=template, title="3D computed jet · drag to rotate; slider selects a saved time",
                         height=700, scene=dict(xaxis_title="x [mm]", yaxis_title="y [mm]", zaxis_title="z [mm]",
                         zaxis=dict(autorange="reversed"), aspectmode="data"), margin=dict(t=60, b=90))
    figures.append(flight)

    if len(trajectory):
        final = trajectory[trajectory["time_s"] == trajectory["time_s"].max()]
        live = final[final["deposited"] < .5]
        if len(live):
            states = make_subplots(rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "scene"}]],
                                  subplot_titles=("Active jet: prescribed convective surface charge", "Active jet: material temperature"))
            for col, values, units, colorscale in (
                (1, live["charge_C"] / np.maximum(live["segment_length_m"] * 2 * np.pi * live["radius_m"], 1e-300) * 1e6, "µC/m²", "RdBu_r"),
                (2, live["temperature_K"] - 273.15, "°C", "Inferno"),
            ):
                states.add_trace(go.Scatter3d(x=live["x_m"] * 1e3, y=live["y_m"] * 1e3, z=live["z_m"] * 1e3,
                                              mode="markers", marker=dict(size=4, color=values, colorscale=colorscale,
                                              colorbar=dict(title=units, x=.45 if col == 1 else 1.0, len=.75)),
                                              name=units, hovertemplate="x=%{x:.4g} mm<br>y=%{y:.4g} mm<br>z=%{z:.4g} mm<br>value=%{marker.color:.4g}<extra></extra>"), row=1, col=col)
            axis = dict(xaxis_title="x [mm]", yaxis_title="y [mm]", zaxis_title="z [mm]", zaxis=dict(autorange="reversed"), aspectmode="data")
            states.update_layout(template=template, height=550, scene=axis, scene2=axis, showlegend=False)
            figures.append(states)

    collector = make_subplots(rows=2, cols=2, subplot_titles=(
        "Resolved fiber centerline chords on collector", "Polymer areal mass [mg/m²]",
        "Impact diameters [µm] — parcel-count weighted", "Chord orientation [deg] — length weighted"))
    segments, masses, lengths = deposited_segments(deposits)
    xs, ys = _line_xy(segments)
    collector.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(width=1.5, color="#4f46e5"), name="Consecutive landed IDs"), row=1, col=1)
    if len(deposits):
        collector.add_trace(go.Scatter(x=deposits["collector_x_m"] * 1e3, y=deposits["collector_y_m"] * 1e3,
                                       mode="markers", marker=dict(size=4, color="#f97316"), name="Impact parcels"), row=1, col=1)
        collector.add_trace(go.Histogram(x=deposits["radius_m"] * 2e6, nbinsx=30, name="Wet at impact", opacity=.65, marker_color="#6366f1"), row=2, col=1)
        collector.add_trace(go.Histogram(x=deposits["dry_radius_m"] * 2e6, nbinsx=30, name="Dry at unchanged length", opacity=.55, marker_color="#f97316"), row=2, col=1)
    collector.add_trace(go.Heatmap(x=(coverage["xedges"][1:] + coverage["xedges"][:-1]) * 500,
                                   y=(coverage["yedges"][1:] + coverage["yedges"][:-1]) * 500,
                                   z=coverage["areal_polymer_mass_kg_m2"] * 1e6, colorscale="Magma", zmin=0,
                                   zmax=max(float(coverage["areal_polymer_mass_kg_m2"].max() * 1e6), 1e-30),
                                   colorbar=dict(title="mg/m²", len=.4, y=.8), name="Mass density"), row=1, col=2)
    if len(segments):
        delta = segments[:, 1] - segments[:, 0]
        angle = (np.degrees(np.arctan2(delta[:, 1], delta[:, 0])) + 90) % 180 - 90
        hist, edges = np.histogram(angle, bins=np.linspace(-90, 90, 37), weights=lengths)
        collector.add_trace(go.Bar(x=(edges[1:] + edges[:-1]) / 2, y=hist * 1e3,
                                   marker_color="#14b8a6", name="Chord length [mm]"), row=2, col=2)
    collector.update_xaxes(title_text="Material x [mm]", row=1, col=1)
    collector.update_yaxes(title_text="Material y [mm]", scaleanchor="x", scaleratio=1, row=1, col=1)
    collector.update_xaxes(title_text="Material x [mm]", row=1, col=2)
    collector.update_yaxes(title_text="Material y [mm]", row=1, col=2)
    collector.update_xaxes(title_text="Diameter [µm]", row=2, col=1)
    collector.update_yaxes(title_text="Parcel count", row=2, col=1)
    collector.update_xaxes(title_text="Orientation modulo 180°", row=2, col=2)
    collector.update_yaxes(title_text="Centerline length [mm]", row=2, col=2)
    collector.update_layout(template=template, height=850, barmode="overlay", legend=dict(orientation="h", y=-.15), margin=dict(b=110))
    figures.append(collector)

    if len(field):
        r, z, phi = field_grid(field, "Phi_V")
        _, _, er = field_grid(field, "E_r_V_m")
        _, _, ez = field_grid(field, "E_z_V_m")
        mag = np.hypot(er, ez)
        fieldfig = make_subplots(rows=2, cols=2, specs=[[{}, {}], [{"type": "surface"}, {"type": "scene"}]],
                                subplot_titles=("Potential [V]", "log10 |E| [V/m]",
                                                "Potential surface — r,z,Φ", "Needle surface charge [µC/m²]"))
        fieldfig.add_trace(go.Heatmap(x=r * 1e3, y=z * 1e3, z=phi, colorscale="RdBu_r",
                                      colorbar=dict(title="V", len=.4, x=.45, y=.8)), row=1, col=1)
        logmag = np.full_like(mag, np.nan)
        np.log10(mag, out=logmag, where=mag > 0)
        fieldfig.add_trace(go.Heatmap(x=r * 1e3, y=z * 1e3, z=logmag, colorscale="Viridis",
                                      colorbar=dict(title="log₁₀ E", len=.4, x=1., y=.8)), row=1, col=2)
        ir = np.unique(np.linspace(0, len(r) - 1, min(96, len(r))).astype(int))
        iz = np.unique(np.linspace(0, len(z) - 1, min(100, len(z))).astype(int))
        fieldfig.add_trace(go.Surface(x=r[ir] * 1e3, y=z[iz] * 1e3, z=phi[np.ix_(iz, ir)] / 1e3,
                                      colorscale="RdBu_r", showscale=False, name="Potential"), row=2, col=1)
        needle = electrode_charge(run_dir)
        if len(needle):
            selected = np.unique(np.linspace(0, len(needle) - 1, min(600, len(needle))).astype(int))
            theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
            faces = needle[selected]
            xx = (faces["r_m"][:, None] * np.cos(theta)).ravel() * 1e3
            yy = (faces["r_m"][:, None] * np.sin(theta)).ravel() * 1e3
            zz = np.repeat(faces["z_m"], len(theta)) * 1e3
            sigma = np.repeat(faces["sigma_C_m2"], len(theta)) * 1e6
            fieldfig.add_trace(go.Scatter3d(x=xx, y=yy, z=zz, mode="markers", marker=dict(size=2.5, color=sigma,
                                           colorscale="RdBu_r", colorbar=dict(title="µC/m²", len=.35, y=.2)),
                                           name="Revolved electrode faces", hovertemplate="x=%{x:.4g} mm<br>y=%{y:.4g} mm<br>z=%{z:.4g} mm<br>σ=%{marker.color:.4g} µC/m²<extra></extra>"), row=2, col=2)
        fieldfig.update_xaxes(title_text="r [mm]", row=1)
        fieldfig.update_yaxes(title_text="z [mm]", autorange="reversed", row=1)
        fieldfig.update_layout(template=template, height=950, scene=dict(xaxis_title="r [mm]", yaxis_title="z [mm]", zaxis_title="Φ [kV]"),
                               scene2=dict(xaxis_title="x [mm]", yaxis_title="y [mm]", zaxis_title="z [mm]", aspectmode="data"), showlegend=False)
        figures.append(fieldfig)

    if len(diagnostics):
        checks = make_subplots(rows=1, cols=2, subplot_titles=("Numerical conservation errors", "Resolved dynamics"))
        for key, label in (("mass_relative_error", "Mass"), ("polymer_relative_error", "Polymer"), ("charge_relative_error", "Charge")):
            checks.add_trace(go.Scatter(x=diagnostics["time_s"], y=np.maximum(np.abs(diagnostics[key]), 1e-18), name=label), row=1, col=1)
        checks.add_trace(go.Scatter(x=diagnostics["time_s"], y=diagnostics["max_speed_m_s"], name="Max speed [m/s]"), row=1, col=2)
        checks.update_xaxes(title_text="t [s]")
        checks.update_yaxes(type="log", title_text="Relative error", row=1, col=1)
        checks.update_yaxes(title_text="m/s", row=1, col=2)
        checks.update_layout(template=template, height=380)
        figures.append(checks)

    plasma_path = run_dir / "field/plasma_final.csv"
    if plasma_path.exists():
        plasma = read_table(plasma_path)
        if len(plasma):
            plasmafig = make_subplots(rows=1, cols=2, subplot_titles=("Gas electron density log10 [m⁻³]", "Mean electron energy [eV]"))
            for col, name in enumerate(("ne_m3", "mean_energy_eV"), 1):
                rr, zz, values = field_grid(plasma, name)
                if name == "ne_m3":
                    values = np.log10(np.maximum(values, 1.))
                plasmafig.add_trace(go.Heatmap(x=rr * 1e3, y=zz * 1e3, z=values, colorscale="Inferno", showscale=False), row=1, col=col)
            plasmafig.update_layout(template=template, height=420, title="Optional gas-plasma diagnostic — separate axisymmetric model")
            plasmafig.update_xaxes(title_text="r [mm]")
            plasmafig.update_yaxes(title_text="z [mm]", autorange="reversed")
            figures.append(plasmafig)

    tiles = [
        ("Landed parcels", str(summary["deposited_parcels"])),
        ("Connected chords", str(summary["connected_deposited_chords"])),
        ("Mapped polymer", f"{summary['mapped_polymer_mass_kg'] * 1e9:.4g} µg"),
        ("Nematic order", "n/a" if summary["orientation_nematic_order_0_to_1"] is None else f"{summary['orientation_nematic_order_0_to_1']:.4f}"),
    ]
    cards = "".join(f"<div class='card'><small>{html.escape(label)}</small><strong>{html.escape(value)}</strong></div>" for label, value in tiles)
    sections = "\n".join(fig.to_html(full_html=False, include_plotlyjs=True if index == 0 else False,
                                   config={"responsive": True, "displaylogo": False, "toImageButtonOptions": {"format": "svg"}})
                          for index, fig in enumerate(figures))
    notes = "".join(f"<li>{html.escape(note)}</li>" for note in summary["notes"])
    warnings = "".join(f"<li>{html.escape(note)}</li>" for note in summary.get("applicability_warnings", []))
    document = """<!doctype html><html lang="ru"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Electrospinning · computational report</title><style>
body{margin:0;background:#f4f6fc;color:#172033;font-family:system-ui,sans-serif}main{max-width:1450px;margin:auto;padding:28px}
h1{font-size:32px;margin-bottom:8px}p,li{line-height:1.6}.sub{color:#58657a}.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:25px 0}
.card{background:white;border:1px solid #e2e8f0;border-radius:14px;padding:18px}.card small{display:block;color:#64748b}.card strong{display:block;font-size:26px;margin-top:8px}
.plots{background:white;border-radius:16px;padding:10px}.note{padding:18px;background:#e8edf8;border-radius:12px;margin:20px 0}details{background:white;padding:18px;border-radius:12px}pre{overflow:auto;font-size:12px}@media(max-width:750px){main{padding:12px}.cards{grid-template-columns:repeat(2,1fr)}}
</style><main><h1>Electrospinning · расчётный отчёт</h1>
<p class="sub">Сохранённые результаты модели: 3D-динамика струи, осаждение на коллекторе и электрическое поле. Все графики работают без интернета.</p>
<div class="note">Карта коллектора показывает рассчитанные отрезки между соседними осевшими узлами. Контакты волокон, укладка по толщине и пористость здесь не разрешаются. Параметры материала требуют экспериментальной калибровки.</div>
""" + f"<div class='cards'>{cards}</div><div class='plots'>{sections}</div><h2>Как читать результаты</h2><ul>{notes}</ul><h2>Границы применимости</h2><ul>{warnings}</ul>" + \
        "<details><summary>Параметры расчёта</summary><pre>" + html.escape(json.dumps(config, ensure_ascii=False, indent=2)) + "</pre></details>" + \
        "<p class='sub'>Численные данные: trajectory.csv, deposition.csv, diagnostics.csv, coverage.csv, analysis.json. Итоговая цепочка для ParaView: trajectory.vtp. Статический рисунок: report.png.</p></main></html>"
    html_path = run_dir / "report.html"
    html_path.write_text(document, encoding="utf-8")
    image_path = _static_report(run_dir, trajectory, deposits, field, coverage)
    paths = [html_path, image_path, run_dir / "analysis.json", run_dir / "coverage.csv"]
    for name in ("trajectory.vtp", "final_jet.csv", "emitter_charge.csv", "field/electrode_charge.csv"):
        if (run_dir / name).exists():
            paths.append(run_dir / name)
    if (run_dir / "field/electrode_charge.csv").exists():
        from .needle_view import render_needle
        needle_paths = render_needle(run_dir)
        if needle_paths:
            needle_link = "<p><a href='needle3d.html'>3D-заряд на стальной игле и жидком мениске: отдельная интерактивная реконструкция</a> · <a href='needle3d.png'>Рисунок PNG</a></p>"
            html_path.write_text(document.replace("<div class='cards'>", needle_link + "<div class='cards'>", 1), encoding="utf-8")
            paths.extend(needle_paths)
    return paths
