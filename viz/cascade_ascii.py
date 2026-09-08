"""
Braille anatomy: project every soma onto the terminal at 2x4 dots per character
and animate the cascade over it.

  view      dorsal: left/right = animal's left/right (from somaSide), top = brain,
            bottom = ventral nerve cord (axis from brain vs VNC soma centroids)
  aspect    corrected for 1:2 character cells, so proportions are anatomical
  layers    dim braille = soma density (the shape of the CNS)
            bright braille (colour by superclass) = neurons spiking in this 5 ms bin
            letters  G giant fiber  T TTMn  P PSI  M MN9  (drawn when they spike)

Run:  python viz/cascade_ascii.py [trial] [data_dir] [--path] [--no-color] [--width N] [--height N]
  --path      only neurons that spike in the window, no background
  --synapses  draw arbors of active neurons (all their synaptic sites, dim) and light
              the synapses onto downstream active neurons when the presynaptic cell
              spikes (bright). Requires viz/fetch_cascade_synapses.py output.
"""
import sys, time, os
import numpy as np
import pandas as pd
from pathlib import Path

argv = sys.argv[1:]
def popflag(name, default):
    if name in argv:
        i = argv.index(name); v = int(argv[i + 1]); del argv[i:i + 2]; return v
    return default
W = popflag("--width", None); H = popflag("--height", None)
args = [a for a in argv if not a.startswith("--")]
flags = [a for a in argv if a.startswith("--")]
TRIAL = int(args[0]) if args and args[0] else None
DATA = args[1] if len(args) > 1 else "data"
PATH_ONLY = "--path" in flags
SYN = "--synapses" in flags
TTY = sys.stdout.isatty()
COLOR = "--no-color" not in flags and TTY
if TTY:
    ts = os.get_terminal_size(); W = W or min(ts.columns - 1, 200); H = H or max(ts.lines - 4, 20)
else:
    W = W or 120; H = H or 48
BIN_MS, WINDOW_MS, FADE = 5, 120, 3
RES = Path("results/gf_escape")

neurons = pd.read_parquet(f"{DATA}/neurons.parquet").dropna(subset=["x", "y", "z"]).copy()
sp = pd.read_parquet(RES / "spikes_real.parquet")
trials = pd.read_csv(RES / "real_electrical.csv")
if TRIAL is None:
    TRIAL = int((trials.gf - 1.5).abs().idxmin())
onset = 300 + TRIAL * 360

# ---- anatomical axes
X = neurons[["x", "y", "z"]].values.astype(float)
sc = neurons["superclass"].fillna("")
brain = X[sc.str.startswith(("cb_", "ol_")).values].mean(0)
vnc = X[sc.str.startswith("vnc").values].mean(0)
ap = vnc - brain; ap /= np.linalg.norm(ap)                       # brain -> VNC, drawn top -> bottom
side = neurons["somaSide"].values
lr = X[side == "R"].mean(0) - X[side == "L"].mean(0)
lr -= ap * (lr @ ap); lr /= np.linalg.norm(lr)                    # mediolateral, orthogonal to ap
u = (X - brain) @ lr        # horizontal (animal's left -> right)
v = (X - brain) @ ap        # vertical (brain -> VNC)

# ---- dot grid with 1:2 cell aspect (braille 2 wide x 4 tall per cell)
DW, DH = W * 2, H * 4
ulo, uhi = np.percentile(u, [0.3, 99.7]); vlo, vhi = np.percentile(v, [0.3, 99.7])
span_u, span_v = uhi - ulo, vhi - vlo
# physical scale: one cell is ~1:2, dots are 2x4 so dots are ~square -> uniform scale in dot units
scale = min(DW / span_u, DH / span_v)
dx = ((u - ulo) * scale).astype(int); dy = ((v - vlo) * scale).astype(int)
ox = (DW - int(span_u * scale)) // 2
dx = np.clip(dx + ox, 0, DW - 1); dy = np.clip(dy, 0, DH - 1)
neurons["dx"] = dx; neurons["dy"] = dy
pos = neurons.set_index("bodyId")[["dx", "dy"]]

def project(xyz):
    uu = (xyz - brain) @ lr; vv = (xyz - brain) @ ap
    px = np.clip(((uu - ulo) * scale).astype(int) + ox, 0, DW - 1)
    py = np.clip(((vv - vlo) * scale).astype(int), 0, DH - 1)
    return px, py

arbor = np.zeros((DH, DW), bool); links = None
if SYN:
    sites = pd.read_parquet(RES / "cascade_sites.parquet")
    ax, ay = project(sites[["x", "y", "z"]].values.astype(float)); arbor[ay, ax] = True
    links = pd.read_parquet(RES / "cascade_links.parquet")
    lx, ly = project(links[["x_pre", "y_pre", "z_pre"]].values.astype(float))
    links = links.assign(dx=lx, dy=ly)
    print(f"arbors: {len(sites):,} synaptic sites; {len(links):,} transmission sites between active neurons")

dens = np.zeros((DH, DW), int); np.add.at(dens, (dy, dx), 1)
thr = np.quantile(dens[dens > 0], 0.35) if (dens > 0).any() else 1
bgdots = dens >= thr

BR = [[0x01, 0x08], [0x02, 0x10], [0x04, 0x20], [0x40, 0x80]]  # braille bit for (row, col) in cell
def cell_char(mask, cy, cx):
    bits = 0
    for r in range(4):
        for c in range(2):
            if mask[cy * 4 + r, cx * 2 + c]:
                bits |= BR[r][c]
    return chr(0x2800 + bits) if bits else " "

PAL = {"ol_intrinsic": 36, "ol_sensory": 36, "visual_projection": 96, "visual_centrifugal": 96,
       "cb_intrinsic": 33, "cb_sensory": 93, "cb_motor": 91, "descending_neuron": 35,
       "ascending_neuron": 95, "vnc_intrinsic": 32, "vnc_motor": 31, "vnc_sensory": 92}
MARK = {"DNp01": "G", "TTMn": "T", "PSI": "P", "MN9": "M"}
DIM = "\033[2;37m"; RESET = "\033[0m"
def col(code, s, bold=False):
    return f"\033[{'1;' if bold else ''}{code}m{s}{RESET}" if COLOR else s

bglayer = arbor if SYN else bgdots
base = [[(DIM + cell_char(bglayer, cy, cx) + RESET if COLOR else cell_char(bglayer, cy, cx))
         if not PATH_ONLY and cell_char(bglayer, cy, cx) != " " else " " for cx in range(W)] for cy in range(H)]
if SYN and not PATH_ONLY:   # faint CNS outline under the arbors
    for cy in range(H):
        for cx in range(W):
            if base[cy][cx] == " ":
                ch = cell_char(bgdots, cy, cx)
                if ch != " ":
                    base[cy][cx] = (DIM + ch + RESET) if COLOR else "."

win = sp[(sp.t_ms >= onset) & (sp.t_ms < onset + WINDOW_MS)].copy()
win = win[win.bodyId.isin(pos.index)]
win["b"] = ((win.t_ms - onset) // BIN_MS).astype(int)
meta = neurons.set_index("bodyId")
print(f"trial {TRIAL}: {len(win)} spikes from {win.bodyId.nunique()} neurons; "
      f"{len(neurons):,} somas at {DW}x{DH} dots ({W}x{H} cells), dorsal view")
time.sleep(1.2)

for b in range(WINDOW_MS // BIN_MS):
    grid = [r[:] for r in base]
    act = np.zeros((DH, DW), bool); actsc = {}
    letters = {}
    for age in range(FADE, -1, -1):
        cur = win[win.b == b - age]
        for bid, t in zip(cur.bodyId, cur.type):
            x, y = pos.at[bid, "dx"], pos.at[bid, "dy"]
            if t in MARK:
                letters[(y // 4, x // 2)] = (MARK[t], meta.at[bid, "superclass"])
            else:
                act[y, x] = True; actsc[(y // 4, x // 2)] = (meta.at[bid, "superclass"], age)
    if links is not None:
        firing = set(win[win.b == b].bodyId)
        if firing:
            L = links[links.bodyId_pre.isin(firing)]
            syn = np.zeros((DH, DW), bool); syn[L.dy.values, L.dx.values] = True
            for cy in range(H):
                for cx in range(W):
                    ch = cell_char(syn, cy, cx)
                    if ch != " ":
                        grid[cy][cx] = col(93, ch, bold=True)      # bright yellow = synapses transmitting now
    for (cy, cx), (s, age) in actsc.items():
        ch = cell_char(act, cy, cx)
        if ch != " ":
            grid[cy][cx] = col(PAL.get(s, 37), ch, bold=(age == 0))
    for (cy, cx), (ch, s) in letters.items():
        grid[cy][cx] = col(PAL.get(s, 37), ch, bold=True)
    if TTY:
        print("\033[H\033[J", end="")
    print(f" t = {b*BIN_MS:>3d} ms after loom onset   dorsal view: brain top, VNC bottom, animal's left on the left"
          f"   [G GF  T TTMn  P PSI  M MN9{'  yellow = synapses transmitting' if SYN else ''}]")
    for r in grid:
        print("".join(r))
    sys.stdout.flush()
    time.sleep(0.25 if TTY else 0)
print("done")
