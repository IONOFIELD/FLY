"""
3D cascade with real morphology.

  skeletons     every neuron active in the window, drawn as its neuprint skeleton (cached in
                data/skeletons/), coloured by superclass; dim until it spikes, bright for a few
                frames after, then dim again
  synapses      transmission sites between active neurons flash yellow when the pre cell spikes
                (from viz/fetch_cascade_synapses.py; skipped if absent)
  outline       neuropil meshes from neuprint if available (cached data/meshes/), else the soma
                cloud, drawn translucent
  controls      play/pause, time slider (5 ms), camera presets dorsal / lateral / anterior

Run:  python viz/cascade_3d.py [trial] [data_dir] [--no-skeletons]
Needs NEUPRINT_TOKEN for the first run (skeletons, meshes); cached afterwards.
Writes results/gf_escape/cascade_3d.html
"""
import os, sys
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pathlib import Path

argv = sys.argv[1:]
flags = [a for a in argv if a.startswith("--")]; args = [a for a in argv if not a.startswith("--")]
TRIAL = int(args[0]) if args and args[0] else None
DATA = args[1] if len(args) > 1 else "data"
RES = Path("results/gf_escape")
BIN_MS, WINDOW_MS, FADE = 5, 120, 4
SK = Path("data/skeletons"); SK.mkdir(parents=True, exist_ok=True)
MESH = Path("data/meshes"); MESH.mkdir(parents=True, exist_ok=True)
PAL = {"ol_intrinsic": "#4dd0e1", "ol_sensory": "#4dd0e1", "visual_projection": "#b39ddb", "visual_centrifugal": "#b39ddb",
       "cb_intrinsic": "#ffd54f", "cb_sensory": "#ffe082", "cb_motor": "#ff8a65", "descending_neuron": "#f06292",
       "ascending_neuron": "#ce93d8", "vnc_intrinsic": "#81c784", "vnc_motor": "#e57373", "vnc_sensory": "#aed581"}

neurons = pd.read_parquet(f"{DATA}/neurons.parquet")
meta = neurons.set_index("bodyId")
sp = pd.read_parquet(RES / "spikes_real.parquet")
trials = pd.read_csv(RES / "real_electrical.csv")
if TRIAL is None:
    TRIAL = int((trials.gf - 1.5).abs().idxmin())
onset = 300 + TRIAL * 360
win = sp[(sp.t_ms >= onset) & (sp.t_ms < onset + WINDOW_MS)].copy()
win["b"] = ((win.t_ms - onset) // BIN_MS).astype(int)
active = sorted(set(win.bodyId) | set(meta.index[meta.type.isin(["DNp01", "TTMn", "PSI"])]))
print(f"trial {TRIAL}: {len(win)} spikes, {len(active)} neurons to draw")

# ---------------------------------------------------------------- skeletons (cached)
client = None
def get_client():
    global client
    if client is None:
        from neuprint import Client
        tok = os.environ.get("NEUPRINT_TOKEN")
        if not tok:
            raise RuntimeError("NEUPRINT_TOKEN not set")
        client = Client("neuprint.janelia.org", dataset="male-cns:v1.0", token=tok)
    return client

def skeleton(bid):
    f = SK / f"{bid}.parquet"
    if f.exists():
        return pd.read_parquet(f)
    from neuprint import fetch_skeleton
    sk = fetch_skeleton(int(bid), format="pandas", client=get_client())
    sk.to_parquet(f, index=False)
    return sk

def skeleton_lines(sk):
    """x,y,z arrays with None separators so one trace draws all branches."""
    P = sk.set_index("rowId")[["x", "y", "z"]]
    seg = sk[sk.link > 0]
    a = P.loc[seg.rowId].values; b = P.loc[seg.link].values
    xs = np.empty(len(a) * 3, object); ys = xs.copy(); zs = xs.copy()
    xs[0::3], xs[1::3], xs[2::3] = a[:, 0], b[:, 0], None
    ys[0::3], ys[1::3], ys[2::3] = a[:, 1], b[:, 1], None
    zs[0::3], zs[1::3], zs[2::3] = a[:, 2], b[:, 2], None
    return xs, ys, zs

skels = {}
if "--no-skeletons" not in flags:
    for k, bid in enumerate(active):
        try:
            skels[bid] = skeleton_lines(skeleton(bid))
        except Exception as e:
            if k == 0:
                print(f"  skeletons unavailable ({e}); falling back to soma markers")
            break
    print(f"  {len(skels)} skeletons")

# ---------------------------------------------------------------- outline: meshes or soma cloud
outline_traces = []
def try_meshes():
    try:
        from neuprint import fetch_roi_mesh
    except ImportError:
        return 0
    rois = ["CentralBrain", "OL(R)", "OL(L)", "VNC", "GNG"]
    got = 0
    for roi in rois:
        f = MESH / f"{roi.replace('(', '_').replace(')', '')}.npz"
        try:
            if f.exists():
                d = np.load(f); v, fc = d["v"], d["f"]
            else:
                m = fetch_roi_mesh(roi, export_path=None, client=get_client())
                v, fc = np.asarray(m.vertices), np.asarray(m.faces)
                np.savez_compressed(f, v=v, f=fc)
            outline_traces.append(go.Mesh3d(x=v[:, 0], y=v[:, 1], z=v[:, 2], i=fc[:, 0], j=fc[:, 1], k=fc[:, 2],
                                            color="#8899aa", opacity=0.08, name=roi, hoverinfo="skip", showlegend=False))
            got += 1
        except Exception:
            continue
    return got
def try_local_meshes():
    got = 0
    for f in sorted(MESH.glob("*.npz")):
        d = np.load(f); v, fc = d["v"], d["f"]
        outline_traces.append(go.Mesh3d(x=v[:, 0], y=v[:, 1], z=v[:, 2], i=fc[:, 0], j=fc[:, 1], k=fc[:, 2],
                                        color="#8899aa", opacity=0.07, name=f.stem, hoverinfo="skip", showlegend=False))
        got += 1
    return got


if not (try_local_meshes() or try_meshes()):
    bg = neurons.dropna(subset=["x", "y", "z"]).sample(min(30000, len(neurons)), random_state=0)
    outline_traces.append(go.Scatter3d(x=bg.x, y=bg.y, z=bg.z, mode="markers", name="all somas",
                                       marker=dict(size=1, color="#8899aa", opacity=0.12), hoverinfo="skip", showlegend=False))
    print("  outline: soma cloud (no ROI meshes available)")
else:
    print(f"  outline: neuropil meshes ({len(outline_traces)})")

# ---------------------------------------------------------------- synapses between active neurons
links = None
if (RES / "cascade_links.parquet").exists():
    links = pd.read_parquet(RES / "cascade_links.parquet")
    print(f"  {len(links)} transmission synapses")

# ---------------------------------------------------------------- traces: one per neuron, colour animated
spike_bins = win.groupby("bodyId").b.apply(set).to_dict()
def brightness(bid, b):
    bins = spike_bins.get(bid, set())
    ages = [b - s for s in bins if 0 <= b - s <= FADE]
    return 1.0 - min(ages) / (FADE + 1) if ages else 0.0

def hex_alpha(hx, a):
    r, g, bl = int(hx[1:3], 16), int(hx[3:5], 16), int(hx[5:7], 16)
    return f"rgba({r},{g},{bl},{a:.2f})"

neuron_traces = []
seen_groups = set()
for bid in active:
    sc = meta.at[bid, "superclass"] or "other"; col = PAL.get(sc, "#cccccc"); t = meta.at[bid, "type"]
    first = sc not in seen_groups; seen_groups.add(sc)
    common = dict(name=sc, legendgroup=sc, showlegend=first, hovertext=f"{t} {bid}", hoverinfo="text")
    if bid in skels:
        xs, ys, zs = skels[bid]
        neuron_traces.append(go.Scatter3d(x=xs, y=ys, z=zs, mode="lines",
                                          line=dict(color=hex_alpha(col, 0.15), width=2), **common))
    else:
        r = meta.loc[bid]
        neuron_traces.append(go.Scatter3d(x=[r.x], y=[r.y], z=[r.z], mode="markers",
                                          marker=dict(size=4, color=hex_alpha(col, 0.15)), **common))
syn_trace = go.Scatter3d(x=[], y=[], z=[], mode="markers", name="synapses transmitting now",
                         marker=dict(size=1.6, color="#ffb74d", opacity=0.7), hoverinfo="skip")

frames = []
for b in range(WINDOW_MS // BIN_MS):
    data = []
    for bid, tr in zip(active, neuron_traces):
        a = 0.12 + 0.88 * brightness(bid, b)
        col = PAL.get(meta.at[bid, "superclass"], "#cccccc")
        if tr.mode == "lines":
            data.append(go.Scatter3d(line=dict(color=hex_alpha(col, a), width=2 + 4 * brightness(bid, b))))
        else:
            data.append(go.Scatter3d(marker=dict(color=hex_alpha(col, a), size=4 + 6 * brightness(bid, b))))
    if links is not None:
        firing = set(win[win.b == b].bodyId)
        L = links[links.bodyId_pre.isin(firing)]
        data.append(go.Scatter3d(x=L.x_pre, y=L.y_pre, z=L.z_pre))
    frames.append(go.Frame(data=data, name=str(b * BIN_MS),
                           traces=list(range(len(outline_traces), len(outline_traces) + len(data))),
                           layout=go.Layout(title=f"t = {b*BIN_MS:>3d} ms after loom onset")))

fig = go.Figure(data=outline_traces + neuron_traces + ([syn_trace] if links is not None else []), frames=frames)

# camera presets from the data's anatomical axes
X = neurons.dropna(subset=["x", "y", "z"])[["x", "y", "z"]].values.astype(float)
scc = neurons.dropna(subset=["x", "y", "z"])["superclass"].fillna("")
ap = X[scc.str.startswith("vnc").values].mean(0) - X[scc.str.startswith(("cb_", "ol_")).values].mean(0); ap /= np.linalg.norm(ap)
side = neurons.dropna(subset=["x", "y", "z"])["somaSide"].values
lr = X[side == "R"].mean(0) - X[side == "L"].mean(0); lr -= ap * (lr @ ap); lr /= np.linalg.norm(lr)
dv = np.cross(lr, ap)
def cam(eye, up):
    return dict(eye=dict(x=float(eye[0]), y=float(eye[1]), z=float(eye[2])), up=dict(x=float(up[0]), y=float(up[1]), z=float(up[2])))
cams = {"dorsal": cam(-2.0 * dv, -ap), "lateral": cam(2.0 * lr, -dv), "anterior": cam(-2.0 * ap, -dv)}

fig.update_layout(
    title=f"Loom-evoked cascade, MaleCNS v1.0 LIF, trial {TRIAL}",
    scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
               aspectmode="data", bgcolor="#0b0e14", camera=cams["dorsal"]),
    paper_bgcolor="#0b0e14", font=dict(color="#e0e0e0"),
    legend=dict(title="cell class (click to hide)", itemsizing="constant", font=dict(size=10), x=0.99, y=0.95,
                bgcolor="rgba(11,14,20,0.6)"),
    margin=dict(l=0, r=0, t=70, b=0),
    updatemenus=[
        dict(type="buttons", direction="right", showactive=True, x=0.5, y=1.06, xanchor="center", yanchor="top",
             bgcolor="#1c2331", bordercolor="#3a4a63", font=dict(color="#e0e0e0"),
             buttons=[dict(label=f"view: {n}", method="relayout", args=[{"scene.camera": c}]) for n, c in cams.items()]),
        dict(type="buttons", direction="right", showactive=False, x=0.01, y=0.0, xanchor="left", yanchor="bottom",
             bgcolor="#1c2331", bordercolor="#3a4a63", font=dict(color="#e0e0e0"), buttons=[
            dict(label="play", method="animate", args=[None, dict(frame=dict(duration=150, redraw=True), fromcurrent=True)]),
            dict(label="pause", method="animate", args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate")])])],
    sliders=[dict(steps=[dict(method="animate", label=f.name, args=[[f.name], dict(mode="immediate", frame=dict(duration=0, redraw=True))])
                        for f in frames], currentvalue=dict(prefix="ms after loom onset: ", font=dict(size=12)),
                  x=0.14, len=0.84, y=0.0, yanchor="bottom", pad=dict(t=10), font=dict(size=9))])
out = RES / "cascade_3d.html"
fig.write_html(out, include_plotlyjs="cdn", auto_play=False)
print("wrote", out)
