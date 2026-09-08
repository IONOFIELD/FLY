"""
Braille dorsal/landscape view of the CNS with the saved GF-escape cascade animated over it.

Run:  python viz/cascade_ascii.py [trial] [data_dir] [--synapses] [--path] [--landscape|--portrait]
                                  [--stretch] [--no-color] [--width N] [--height N]
  --synapses  arbors of active neurons (dim) + transmission sites lighting when the pre cell
              spikes (yellow). Requires viz/fetch_cascade_synapses.py output.
  --path      no CNS outline, only active neurons (and arbors if --synapses)
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns.anatomy import Projector

argv = sys.argv[1:]
def popint(name):
    if name in argv:
        i = argv.index(name); v = int(argv[i + 1]); del argv[i:i + 2]; return v
W, H = popint("--width"), popint("--height")
flags = [a for a in argv if a.startswith("--")]; args = [a for a in argv if not a.startswith("--")]
TRIAL = int(args[0]) if args and args[0] else None
DATA = args[1] if len(args) > 1 else "data"
RES = Path("results/gf_escape")

neurons = pd.read_parquet(f"{DATA}/neurons.parquet")
sp = pd.read_parquet(RES / "spikes_real.parquet")
trials = pd.read_csv(RES / "real_electrical.csv")
if TRIAL is None:
    TRIAL = int((trials.gf - 1.5).abs().idxmin())
onset = 300 + TRIAL * 360
land = True if "--landscape" in flags else False if "--portrait" in flags else None
proj = Projector(neurons, W, H, color=("--no-color" not in flags), landscape=land, stretch=("--stretch" in flags))

base, links = None, None
if "--synapses" in flags:
    sites = pd.read_parquet(RES / "cascade_sites.parquet")
    ax, ay = proj.project(sites[["x", "y", "z"]].values.astype(float))
    arbor = np.zeros((proj.DH, proj.DW), bool); arbor[ay, ax] = True
    base = proj.base_grid(arbor, outline=("--path" not in flags))
    links = pd.read_parquet(RES / "cascade_links.parquet")
    lx, ly = proj.project(links[["x_pre", "y_pre", "z_pre"]].values.astype(float))
    links = links.assign(dx=lx, dy=ly)
    print(f"arbors: {len(sites):,} sites; {len(links):,} transmission synapses between active neurons")
elif "--path" in flags:
    base = [[" "] * proj.W for _ in range(proj.H)]
w = sp[(sp.t_ms >= onset) & (sp.t_ms < onset + 120)]
print(f"trial {TRIAL}: {len(w)} spikes from {w.bodyId.nunique()} neurons; {proj.W}x{proj.H} cells, {proj.view_label}")
proj.animate(sp, onset, window_ms=120, title=f"loom trial {TRIAL}", base=base, links=links)
print("done")
