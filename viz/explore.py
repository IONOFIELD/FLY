"""
Interactive: pick a sensory group, drive it, watch the cascade.

Builds ONE network with every selectable group as a Poisson source (about a
minute), then loops: choose group and rate, run 300 ms (stim 100 ms), animate
the braille cascade, print which types and regions fired. Each run ~10 s.

Run:  python viz/explore.py [data_dir] [--sez GAIN] [--landscape|--portrait] [--stretch]
      landscape (brain left, VNC right) is the default on wide terminals
"""
import sys, time
import numpy as np
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns import load_graph, CNSModel, LIFParams
from flycns.anatomy import Projector
from flycns import ascii as A

argv = sys.argv[1:]
SEZ = 1.0
if "--sez" in argv:
    i = argv.index("--sez"); SEZ = float(argv[i + 1]); del argv[i:i + 2]
LANDSCAPE = None
if "--landscape" in argv: LANDSCAPE = True; argv.remove("--landscape")
if "--portrait" in argv: LANDSCAPE = False; argv.remove("--portrait")
STRETCH = "--stretch" in argv
if STRETCH: argv.remove("--stretch")
DATA = argv[0] if argv else "data"

neurons, edges = load_graph(DATA)
t = neurons["type"].fillna("")
sub = neurons["subclass"].fillna("") if "subclass" in neurons else pd.Series("", index=neurons.index)

def types_where(mask):
    return sorted(neurons.loc[mask, "type"].dropna().unique())

GROUPS = {
    "loom (LC4 + LPLC2)":                 types_where(t.isin(["LC4", "LPLC2"])),
    "all loom-sensitive LC types":        types_where(t.isin(["LC4", "LPLC2", "LC6", "LC16", "LC26", "LPLC1", "LC9", "LC17", "LC12"])),
    "photoreceptors R1-R6 (histamine)":   types_where(t == "R1-R6"),
    "lamina L1 + L2 (light OFF proxy)":   types_where(t.isin(["L1", "L2"])),
    "sound: JO-A + JO-B":                 types_where(t.str.match(r"^JO-[AB]") & ~t.str.endswith("unclear")),
    "wind / gravity: JO-C/E/F":           types_where(t.str.match(r"^JO-[CEF]")),
    "pharyngeal gustatory":               types_where(sub == "pharyngeal sensillum") or types_where(t.str.startswith("aPhM")),
    "labellar bristle gustatory":         types_where(sub == "labellar bristle") or types_where(t.str.match(r"^LB\d")),
    "taste peg":                          types_where(sub == "taste peg") or types_where(t.str.endswith("tpGRN")),
    "labellar afferents (BM_Taste)":      types_where(t == "BM_Taste"),
    "leg mechanosensory (claw/ bristle)": types_where(t.str.match(r"^(claw|BM_InOm|LB1)")),
}
GROUPS = {k: v for k, v in GROUPS.items() if v}
stim_all = sorted({x for v in GROUPS.values() for x in v})
print(f"building network with {len(stim_all)} selectable stimulus types (SEZ gain {SEZ})...")
m = CNSModel(neurons, edges, stim_all, LIFParams(region_gains={"SEZ": SEZ, "other": 1.0}), electrical=True)
m.store()
proj = Projector(neurons, landscape=LANDSCAPE, stretch=STRETCH)
print("ready.\n")

while True:
    keys = list(GROUPS)
    for i, k in enumerate(keys):
        print(f"  {i+1:>2}. {k:<38} {len(GROUPS[k]):>3} types, {int(neurons.type.isin(GROUPS[k]).sum()):>5} neurons")
    print("   q. quit")
    choice = input("group> ").strip()
    if choice.lower() in ("q", "quit", ""):
        break
    try:
        grp = keys[int(choice) - 1]
    except (ValueError, IndexError):
        print("pick a number"); continue
    rate = input("rate Hz [50]> ").strip(); rate = float(rate) if rate else 50.0
    types = GROUPS[grp]
    m.restore()
    rates = m.stim["type"].isin(types).values * rate
    m.run(100); on = m.t_ms; m.set_stim_rates(rates); m.run(100); m.set_stim_rates(np.zeros(len(m.stim))); m.run(200)
    sp = m.spike_frame()
    w = sp[(sp.t_ms >= on) & (sp.t_ms < on + 300)]
    print(f"\n{grp} at {rate:.0f} Hz: {len(w)} spikes from {w.bodyId.nunique()} neurons in 300 ms; "
          f"whole-CNS {len(w)/len(m.lif_ids)/0.3:.4f} Hz/neuron")
    input("enter to play cascade...")
    proj.animate(sp, on, window_ms=200, title=f"{grp} @ {rate:.0f} Hz")
    print("\ntop responding types (spikes):")
    print(w.groupby("type").size().sort_values(ascending=False).head(12).to_string())
    A.region_bar(m, [on], window_ms=300)
    for r in ["DNp01", "TTMn", "MN9", "MN11D", "DLMn a,b", "DVMn 1a-c"]:
        ids = set(m.lif_ids[m.readout_index(r)])
        if ids:
            print(f"  {r:<10} {int(w.bodyId.isin(ids).sum())} spikes")
    print()
