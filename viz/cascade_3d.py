"""
3D animation of a loom-evoked cascade through the MaleCNS.

Reads results/gf_escape/spikes_real.parquet and data/neurons.parquet.
Shows soma positions; neurons light up in the time bin they spike, coloured
by superclass, over a grey background of all somas. One frame per BIN_MS.
Writes results/gf_escape/cascade_3d.html (open in a browser).

Usage: python viz/cascade_3d.py [trial_index] [data_dir]
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go

TRIAL = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] else None
DATA = sys.argv[2] if len(sys.argv) > 2 else "data"
RES = Path("results/gf_escape")
BIN_MS, WINDOW_MS, FADE_BINS = 5, 120, 3
BG_N = 25000

neurons = pd.read_parquet(f"{DATA}/neurons.parquet").dropna(subset=["x", "y", "z"])
sp = pd.read_parquet(RES / "spikes_real.parquet")
trials = pd.read_csv(RES / "real_electrical.csv")

# pick the trial with GF count closest to 1.5/cell unless given
if TRIAL is None:
    TRIAL = int((trials.gf - 1.5).abs().idxmin())
# loom onsets are reconstructed from the protocol: gap 300 + k*(300+60)
onset = 300 + TRIAL * 360
print(f"trial {TRIAL}, onset {onset} ms, GF {trials.gf[TRIAL]:.1f}/cell, gain {trials.gain[TRIAL]:.2f}")

win = sp[(sp.t_ms >= onset) & (sp.t_ms < onset + WINDOW_MS)].copy()
win["bin"] = ((win.t_ms - onset) // BIN_MS).astype(int)
win = win.drop(columns=["superclass"]).merge(neurons[["bodyId", "x", "y", "z", "superclass"]], on="bodyId", how="inner")
print(f"{len(win):,} spikes from {win.bodyId.nunique():,} neurons in window")

classes = sorted(win.superclass.dropna().unique())
palette = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0",
           "#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8"]
color = {c: palette[k % len(palette)] for k, c in enumerate(classes)}

rng = np.random.default_rng(0)
bg = neurons.sample(min(BG_N, len(neurons)), random_state=0)
base = go.Scatter3d(x=bg.x, y=bg.y, z=bg.z, mode="markers", name="all somas",
                    marker=dict(size=1.2, color="#bbbbbb", opacity=0.15), hoverinfo="skip")

n_bins = WINDOW_MS // BIN_MS
frames = []
for b in range(n_bins):
    traces = []
    for c in classes:
        recent = win[(win.superclass == c) & (win.bin <= b) & (win.bin > b - FADE_BINS)]
        age = (b - recent.bin).values
        traces.append(go.Scatter3d(
            x=recent.x, y=recent.y, z=recent.z, mode="markers", name=c,
            marker=dict(size=np.clip(6 - 1.5 * age, 2, 6), color=color[c],
                        opacity=0.9), text=recent.type, hoverinfo="text+name"))
    frames.append(go.Frame(data=[base] + traces, name=str(b * BIN_MS),
                           layout=go.Layout(title=f"t = {b*BIN_MS:>3d} ms after loom onset")))

fig = go.Figure(data=frames[0].data, frames=frames)
fig.update_layout(
    title="Loom-evoked cascade, MaleCNS v1.0 LIF (trial %d)" % TRIAL,
    scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
               aspectmode="data", bgcolor="black"),
    paper_bgcolor="black", font=dict(color="white"), legend=dict(itemsizing="constant"),
    updatemenus=[dict(type="buttons", showactive=False, x=0.05, y=0.05, buttons=[
        dict(label="Play", method="animate",
             args=[None, dict(frame=dict(duration=120, redraw=True), fromcurrent=True)]),
        dict(label="Pause", method="animate",
             args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate")])])],
    sliders=[dict(steps=[dict(method="animate", label=f.name,
                              args=[[f.name], dict(mode="immediate", frame=dict(duration=0, redraw=True))])
                        for f in frames], currentvalue=dict(prefix="ms: "))])
out = RES / "cascade_3d.html"
fig.write_html(out, include_plotlyjs="cdn", auto_play=False)
print("wrote", out)
