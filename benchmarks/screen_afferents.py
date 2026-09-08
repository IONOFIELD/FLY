"""
Afferent sufficiency screen for the feeding circuit (after Shiu et al. 2024,
who screened 106 SEZ cell types for MN9 activation at 50 Hz).

Every sensory-superclass type with >= MIN_N neurons and any path within two
hops of MN9 is driven alone at 50 Hz for 200 ms x N_TRIALS. We record MN9
spikes per cell, whole-CNS rate, and MN9 laterality. Types that drive MN9 are
sugar-like candidates; a second pass co-applies each remaining type with the
strongest driver to find suppressors (bitter-like). Results are a table, not
an annotation: modality is inferred from wiring + dynamics only.

Run:  python benchmarks/screen_afferents.py [n_trials] [data_dir] [hz] [w_scale]
      e.g. 3 data 200 0.3   (rate test)     3 data 50 1.0   (Shiu gain test)
      python benchmarks/screen_afferents.py 3 data 50 0.5 aPhM2a,aPhM5,TPMN1   (subset)
Writes results/feeding/afferent_screen.csv and prints suggested drive sets.
"""
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns import load_graph, CNSModel, LIFParams
from flycns.bench import pulse_protocol, readout_rates

N_TRIALS = int(sys.argv[1]) if len(sys.argv) > 1 else 3
DATA = sys.argv[2] if len(sys.argv) > 2 else "data"
HZ = float(sys.argv[3]) if len(sys.argv) > 3 else 50          # drive rate
W_SCALE = float(sys.argv[4]) if len(sys.argv) > 4 else 0.3    # diagnostic only; suite uses 0.3
ONLY = sys.argv[5].split(",") if len(sys.argv) > 5 else None     # restrict to these types
MIN_N, WINDOW = 4, 250
OUT = Path("results/feeding"); OUT.mkdir(parents=True, exist_ok=True)

neurons, edges = load_graph(DATA)
meta = neurons.set_index("bodyId")
infl = pd.read_csv(OUT / "afferent_influence.csv", index_col=0) if (OUT / "afferent_influence.csv").exists() else None
if infl is None:
    raise SystemExit("run benchmarks/diagnose_feeding.py first")
cands = [t for t in infl.index if infl.loc[t, "n_type"] >= MIN_N]
if ONLY:
    cands = [t for t in ONLY if t in infl.index]
print(f"screening {len(cands)} afferent types at {HZ} Hz, w_scale {W_SCALE}, {N_TRIALS} trials each")
tag = f"_{int(HZ)}Hz_w{W_SCALE}"


def drive(types, rates):
    m = CNSModel(neurons, edges, types, LIFParams(w_scale=W_SCALE), electrical=True)
    on = pulse_protocol(m, rates, n_trials=N_TRIALS, pulse_ms=200)
    r = readout_rates(m, on, ["MN9"], WINDOW, meta)
    return dict(mn9=r.spikes_per_cell.mean() if len(r) else 0.0,
                mn9_L=r[r.side == "L"].spikes_per_cell.mean() if (r.side == "L").any() else np.nan,
                mn9_R=r[r.side == "R"].spikes_per_cell.mean() if (r.side == "R").any() else np.nan,
                cns=m.population_rate_hz())


rows = []
t0 = time.time()
for k, t in enumerate(cands):
    s = drive([t], {t: HZ})
    rows.append(dict(type=t, n=int(infl.loc[t, "n_type"]), **s))
    print(f"[{k+1:>2}/{len(cands)}] {t:<18} n={rows[-1]['n']:<4} MN9 {s['mn9']:6.2f}  CNS {s['cns']:.4f} Hz  ({time.time()-t0:.0f}s)")
scr = pd.DataFrame(rows).sort_values("mn9", ascending=False)

# suppression pass: co-apply each non-driver with the strongest driver
drivers = scr[(scr.mn9 > 0.5) & (scr.cns < 0.05)]
if len(drivers):
    top = drivers.iloc[0]["type"]
    base = float(drivers.iloc[0]["mn9"])
    supp = []
    for t in scr[scr.mn9 <= 0.5].type:
        s = drive([top, t], {top: HZ, t: HZ})
        supp.append(dict(type=t, mn9_with_driver=s["mn9"], suppression=(base - s["mn9"]) / base))
        print(f"  suppression {t:<18} MN9 {s['mn9']:6.2f} vs {base:.2f} alone  ({(base-s['mn9'])/base:+.2f})")
    scr = scr.merge(pd.DataFrame(supp), on="type", how="left")
scr.to_csv(OUT / f"afferent_screen{tag}.csv", index=False)

print("\nDRIVERS (MN9 > 0.5 spikes/cell, CNS quiet):")
print(drivers[["type", "n", "mn9", "mn9_L", "mn9_R", "cns"]].to_string(index=False))
if "suppression" in scr:
    sup = scr[scr.suppression > 0.3].sort_values("suppression", ascending=False)
    print("\nSUPPRESSORS (>30% reduction of driver response):")
    print(sup[["type", "n", "suppression"]].to_string(index=False))
    print(f"\nsuggested: python benchmarks/feeding.py 6 data {','.join(drivers.type)} {','.join(sup.type)}")
