"""
Draw one neuron's skeleton in braille, in the same anatomical frame as the CNS views.

  lines    skeleton branches (dim); soma 'O'
  --synapses   input sites (cyan) and output sites (yellow) from neuprint
  --partners   also draw the synapses onto the neuron's top downstream partner
  view: dorsal (default), --side (lateral), --front

Run:  python viz/neuron_ascii.py DNp01            (type; picks the left cell, or --right)
      python viz/neuron_ascii.py 10001            (bodyId)
      python viz/neuron_ascii.py TTMn --synapses --side
Needs NEUPRINT_TOKEN. Skeletons are cached in data/skeletons/.
"""
import os, sys
import numpy as np
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns.anatomy import Projector, BR, DIM, RESET

argv = sys.argv[1:]
flags = [a for a in argv if a.startswith("--")]; args = [a for a in argv if not a.startswith("--")]
if not args:
    raise SystemExit(__doc__)
neurons = pd.read_parquet("data/neurons.parquet")
meta = neurons.set_index("bodyId")
if args[0].isdigit():
    bid = int(args[0])
else:
    cands = meta[meta.type == args[0]]
    if not len(cands):
        raise SystemExit(f"no neuron of type {args[0]}")
    side = "R" if "--right" in flags else "L"
    pick = cands[cands.somaSide == side]
    bid = (pick if len(pick) else cands).index[0]
row = meta.loc[bid]
print(f"{row['type']} bodyId {bid} side {row['somaSide']} superclass {row['superclass']} NT {row.get('consensusNt')}")

# ---- skeleton (cached)
cache = Path("data/skeletons"); cache.mkdir(parents=True, exist_ok=True)
f = cache / f"{bid}.parquet"
if f.exists():
    sk = pd.read_parquet(f)
else:
    from neuprint import Client, fetch_skeleton
    c = Client("neuprint.janelia.org", dataset="male-cns:v1.0", token=os.environ["NEUPRINT_TOKEN"])
    sk = fetch_skeleton(bid, format="pandas", client=c)
    sk.to_parquet(f, index=False)
print(f"skeleton: {len(sk):,} nodes")

# ---- densify segments so branches draw as continuous lines
P = sk.set_index("rowId")[["x", "y", "z"]]
seg = sk[sk.link > 0]
a = P.loc[seg.rowId].values.astype(float); b = P.loc[seg.link].values.astype(float)
n_per = np.clip((np.linalg.norm(a - b, axis=1) / 300).astype(int), 1, 40)   # ~one point per 300 nm
pts = np.concatenate([a[i] + (b[i] - a[i]) * np.linspace(0, 1, n_per[i] + 1)[:, None] for i in range(len(a))])

# ---- projection: reuse CNS anatomical axes, but zoom to this neuron
tty = sys.stdout.isatty()
proj = Projector(neurons, landscape=False)      # for axes only
view = "side" if "--side" in flags else "front" if "--front" in flags else "dorsal"
ax_h, ax_v = {"dorsal": (proj.lr, proj.ap), "side": (proj.ap, np.cross(proj.lr, proj.ap)),
              "front": (proj.lr, np.cross(proj.lr, proj.ap))}[view]
def to2d(xyz):
    return (xyz - proj.brain) @ ax_h, (xyz - proj.brain) @ ax_v
u, v = to2d(pts)
if tty:
    ts = os.get_terminal_size(); W, H = ts.columns - 2, ts.lines - 5
else:
    W, H = 110, 45
DW, DH = W * 2, H * 4
pad = 0.03
ulo, uhi = u.min(), u.max(); vlo, vhi = v.min(), v.max()
ulo -= (uhi - ulo) * pad; uhi += (uhi - ulo) * pad; vlo -= (vhi - vlo) * pad; vhi += (vhi - vlo) * pad
scale = min(DW / (uhi - ulo), DH / (vhi - vlo))
ox = (DW - int((uhi - ulo) * scale)) // 2; oy = (DH - int((vhi - vlo) * scale)) // 2
def grid_xy(xyz):
    uu, vv = to2d(xyz)
    return (np.clip(((uu - ulo) * scale).astype(int) + ox, 0, DW - 1),
            np.clip(((vv - vlo) * scale).astype(int) + oy, 0, DH - 1))
px, py = grid_xy(pts)
skel = np.zeros((DH, DW), bool); skel[py, px] = True

layers = [(skel, DIM if tty else "", False)]
if "--synapses" in flags:
    from neuprint import Client, fetch_synapses, NeuronCriteria as NC
    c = Client("neuprint.janelia.org", dataset="male-cns:v1.0", token=os.environ["NEUPRINT_TOKEN"])
    syn = fetch_synapses(NC(bodyId=bid), client=c)
    for kind, colr in [("post", "\033[1;36m"), ("pre", "\033[1;33m")]:
        s = syn[syn.type == kind]
        m = np.zeros((DH, DW), bool)
        if len(s):
            sx, sy = grid_xy(s[["x", "y", "z"]].values.astype(float)); m[sy, sx] = True
        layers.append((m, colr if tty else "", True))
    print(f"synapses: {int((syn.type=='post').sum()):,} inputs (cyan), {int((syn.type=='pre').sum()):,} outputs (yellow)")
if "--partners" in flags:
    from neuprint import Client, fetch_synapse_connections, NeuronCriteria as NC
    c = Client("neuprint.janelia.org", dataset="male-cns:v1.0", token=os.environ["NEUPRINT_TOKEN"])
    edges = pd.read_parquet("data/edges.parquet")
    top = edges[edges.pre == bid].sort_values("weight", ascending=False).iloc[0]
    con = fetch_synapse_connections(NC(bodyId=bid), NC(bodyId=int(top.post)), client=c)
    m = np.zeros((DH, DW), bool)
    sx, sy = grid_xy(con[["x_pre", "y_pre", "z_pre"]].values.astype(float)); m[sy, sx] = True
    layers.append((m, "\033[1;35m" if tty else "", True))
    print(f"magenta: {len(con)} synapses onto top partner {meta.at[int(top.post), 'type']} ({int(top.post)})")

def cell(mask, cy, cx):
    bits = 0
    for r in range(4):
        for cc in range(2):
            if mask[cy * 4 + r, cx * 2 + cc]:
                bits |= BR[r][cc]
    return chr(0x2800 + bits) if bits else " "

soma = grid_xy(np.array([[row.x, row.y, row.z]], dtype=float)) if pd.notna(row.x) else None
print(f"view: {view}   ({W}x{H} cells, {DW}x{DH} dots)   O = soma")
for cy in range(H):
    line = []
    for cx in range(W):
        ch = " "
        for mask, colr, bright in layers:
            c_ = cell(mask, cy, cx)
            if c_ != " ":
                ch = f"{colr}{c_}{RESET}" if colr else c_
        if soma is not None and soma[1][0] // 4 == cy and soma[0][0] // 2 == cx:
            ch = ("\033[1;31mO" + RESET) if tty else "O"
        line.append(ch)
    print("".join(line))
