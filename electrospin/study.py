"""Parameter-study visualization; failures and absent deposition remain explicit."""
from pathlib import Path
import json
import numpy as np


def render_study(out: Path, rows: list[dict], parameters: list[str]) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plotly.graph_objects as go
    # Scalar outcomes vary across studies; use only measured outputs available here.
    candidates = ['deposited_parcels','mean_dry_diameter_same_length_at_hits_m',
                  'delivered_polymer_mass_kg','rms_impact_radius_material_m']
    metrics = [key for key in candidates if any(isinstance(r.get(key),(float,int)) for r in rows)]
    if not metrics:
        keys = set().union(*(r.keys() for r in rows))
        metrics = [k for k in sorted(keys) if k not in parameters and k != 'case' and
                   any(isinstance(r.get(k),(int,float)) and not isinstance(r.get(k),bool) for r in rows)][:4]
    fig,axes = plt.subplots(1,max(1,len(metrics)),figsize=(max(6,4*len(metrics)),4),squeeze=False)
    interactive = go.Figure()
    for j,metric in enumerate(metrics):
        valid = [r for r in rows if r.get('status')=='COMPLETED' and isinstance(r.get(metric),(float,int))]
        xkey = parameters[0]
        x = [r[xkey] for r in valid];y = [r[metric] for r in valid]
        axes[0,j].scatter(x,y,c='#157f91',s=30)
        axes[0,j].set(xlabel=xkey,ylabel=metric);axes[0,j].grid(alpha=.2)
        interactive.add_trace(go.Scatter(x=x,y=y,mode='markers',name=metric,
                             text=[json.dumps({k:r[k] for k in parameters}) for r in valid],
                             hovertemplate='%{text}<br>value=%{y}<extra>%{fullData.name}</extra>'))
    if not metrics:
        axes[0,0].text(.5,.5,'No completed scalar outcomes\nSee study.csv for failure details',ha='center',va='center')
    fig.suptitle('Parameter study — reduced model, illustrative inputs')
    fig.tight_layout();fig.savefig(out/'study.png',dpi=160);plt.close(fig)
    interactive.update_layout(template='plotly_white',title='Parameter study (legend selects outcome; units in column names)',
                              xaxis_title=parameters[0],yaxis_title='Metric value (different units; select one trace)')
    for trace in list(interactive.data)[1:]: trace.visible='legendonly'
    if len(parameters) >= 2 and metrics:
        metric = metrics[1] if len(metrics) > 1 else metrics[0]
        valid = [r for r in rows if r.get('status')=='COMPLETED' and isinstance(r.get(metric),(float,int))]
        spatial = go.Figure(go.Scatter3d(x=[r[parameters[0]] for r in valid],
            y=[r[parameters[1]] for r in valid],z=[r[metric] for r in valid],mode='markers',
            marker=dict(size=6,color=[r[metric] for r in valid],colorscale='Viridis',showscale=True)))
        spatial.update_layout(template='plotly_white',title='Computed response samples (no fitted surface)',
            scene=dict(xaxis_title=parameters[0],yaxis_title=parameters[1],zaxis_title=metric))
        spatial.write_html(out/'response3d.html',include_plotlyjs=True,auto_open=False)
    interactive.write_html(out/'study.html',include_plotlyjs=True,full_html=True,auto_open=False)
