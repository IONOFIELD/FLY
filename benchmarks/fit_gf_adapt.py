"""
Fit the GF spike-triggered adaptation increment against published constraints.

von Reyn et al. 2014: GF produces at most one or two spikes per loom regardless of
stimulus strength, and looms are frequently subthreshold. Constraints scored per
candidate increment b (mV), using the standard loom protocol with lognormal gain:
  C1  spikes per responding GF in [1, 2]           across ALL trials
  C2  spikes per responding GF <= 2 on the strongest quartile of gains (all-or-none)
  C3  response probability in [0.5, 1.0]
  C4  TTMn spikes / GF spikes in [0.8, 1.2]        (relay unaffected)
Rule: choose the SMALLEST b satisfying C1-C4 (minimal modification of the raw model), with at
least 40 trials so the strongest-gain quartile (C2, the binding constraint) is estimable.
Run:  python benchmarks/fit_gf_adapt.py [n_trials] [data_dir]
Writes results/fit_gf_adapt/grid.csv and best.json
"""
import json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns import load_graph, CNSModel, LIFParams, loom_protocol
from flycns.protocols import LOOM_TUNING
from flycns import graph as G

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
if N < 40:
    print(f"WARNING: {N} trials leaves ~{N//4} trials in the strongest-gain quartile, which is the "
          f"binding constraint (C2). Use >= 40.")
DATA = sys.argv[2] if len(sys.argv) > 2 else "data"
OUT = Path("results/fit_gf_adapt"); OUT.mkdir(parents=True, exist_ok=True)
GRID = [0, 5, 10, 15, 20, 30, 45, 60]

neurons, edges = load_graph(DATA)
rows = []
for b in GRID:
    t0 = time.time()
    # override the GF entry of INTRINSIC_OVERRIDES for this candidate
    saved = list(G.INTRINSIC_OVERRIDES)
    G.INTRINSIC_OVERRIDES[:] = [(t, p, (b if t == "DNp01" else v), s) for t, p, v, s in saved]
    import flycns.model as M; M.INTRINSIC_OVERRIDES = G.INTRINSIC_OVERRIDES
    m = CNSModel(neurons, edges, list(LOOM_TUNING), LIFParams(), electrical=True)
    trials = loom_protocol(m, n_trials=N)
    G.INTRINSIC_OVERRIDES[:] = saved
    sp = m.spike_frame()
    gf = set(m.lif_ids[m.readout_index("DNp01")]); tt = set(m.lif_ids[m.readout_index("TTMn")])
    per, hit, gains, ngf, ntt = [], [], [], 0, 0
    for on, gain in trials:
        w = sp[(sp.t_ms >= on) & (sp.t_ms < on + 160)]
        g = w[w.bodyId.isin(gf)]; t = w[w.bodyId.isin(tt)]
        ngf += len(g); ntt += len(t); hit.append(int(len(g) > 0)); gains.append(gain)
        per.append(len(g) / g.bodyId.nunique() if len(g) else np.nan)
    per = np.array(per); gains = np.array(gains); hit = np.array(hit)
    strong = gains >= np.quantile(gains, 0.75)
    r = dict(b_mV=b, per_response=np.nanmean(per), per_response_strong=np.nanmean(per[strong]),
             max_per_response=np.nanmax(per) if np.isfinite(per).any() else 0, hit=hit.mean(),
             ttmn_ratio=(ntt / ngf) if ngf else np.nan, wall_s=round(time.time() - t0))
    r["C1"] = 1.0 <= r["per_response"] <= 2.0
    r["C2"] = r["per_response_strong"] <= 2.0
    r["C3"] = 0.5 <= r["hit"] <= 1.0
    r["C4"] = 0.8 <= (r["ttmn_ratio"] if np.isfinite(r["ttmn_ratio"]) else 0) <= 1.2
    r["ok"] = all(r[c] for c in ["C1", "C2", "C3", "C4"])
    rows.append(r)
    print(f"b={b:>3} mV  per response {r['per_response']:.2f} (strong gains {r['per_response_strong']:.2f}, max {r['max_per_response']:.0f})  "
          f"hit {r['hit']:.2f}  TTMn/GF {r['ttmn_ratio']:.2f}  [{' '.join(c for c in ['C1','C2','C3','C4'] if r[c]):<11}] {'OK' if r['ok'] else ''}  {r['wall_s']}s")
    pd.DataFrame(rows).to_csv(OUT / "grid.csv", index=False)
df = pd.DataFrame(rows); ok = df[df.ok]
best = int(ok.b_mV.min()) if len(ok) else None
print(f"\nsmallest increment satisfying C1-C4: {best} mV" if best is not None else "\nno candidate satisfies all constraints")
(OUT / "best.json").write_text(json.dumps(dict(b_adapt_mV=best, rule="smallest b with C1-C4",
                                               constraints="von Reyn 2014: 1-2 spikes/loom, probabilistic; Tanouye & Wyman 1980 relay"), indent=2))
