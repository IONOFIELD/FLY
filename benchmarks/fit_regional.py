"""
Regional synaptic gain fit (SEZ intrinsic synapses).

The class-wise gain fit (fit_gains.py) showed taste -> MN9 never activates
even at 3x sensory gain, while the local gain feeding needs (~3x) breaks the
escape circuit. Here only synapses from SEZ intrinsic types (GNG, SAD, PRW,
FLA, CAN prefixes) are scaled; everything else stays at unity.
Gustatory afferents are selected by MaleCNS `subclass` (labellar bristle,
taste peg, pharyngeal sensillum) after fetch_annotations.py.
The SEZ population is defined anatomically (>50% of synapses in SEZ compartments,
excluding sensory/motor/efferent) when data/roi_membership.parquet exists; the model
prints which definition it used.

One network per candidate gain set (sensory, relay, local) on top of the shared
w_scale = 0.3. Stimulus sets are swapped with store/restore, so each candidate
costs one build (~50 s) plus a few short protocols. Each candidate is scored on
quick versions of all three benchmarks plus a specificity control:

  E1 GF hit >= 0.8 and 1-2 spikes/responding GF (loom)        [von Reyn 2014]
  E2 TTMn/GF ratio 0.8-1.2                                     [Tanouye & Wyman 1980]
  E3 whole-CNS rate < 0.05 Hz/neuron during loom
  A1 JO alone: GF hit <= 0.3                                    [Pezier & Blagburn 2013]
  A2 JO + loom: GF hit or count > loom alone, or latency shorter
  F1 taste afferents at 50 Hz drive MN9 (> 0.5 spikes/cell)    [Gordon & Scott 2009]
  F2 whole-CNS rate < 0.05 Hz/neuron during taste drive
  S1 JO-FV (wind/auditory) at 50 Hz does NOT drive MN9 (< 0.2) [specificity]

Score = number of checks passed (max 8). Ties broken by smaller total gain
deviation from 1 (prefer the least modified model).
Run:  python benchmarks/fit_regional.py [n_trials] [data_dir]
Writes results/fit_regional/grid.csv and best.json
"""
import json, sys, time, itertools
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns import load_graph, CNSModel, LIFParams
from flycns.protocols import LOOM_TUNING, loom_protocol
from flycns.bench import pulse_protocol, readout_rates, find_types
from flycns import ascii as A

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
DATA = sys.argv[2] if len(sys.argv) > 2 else "data"
OUT = Path("results/fit_regional"); OUT.mkdir(parents=True, exist_ok=True)

SEZ_GRID = [1.5, 1.75, 2.0, 2.25, 2.5]
from flycns.graph import gustatory_afferents
TASTE_FALLBACK = ["aPhM2a", "aPhM5", "PhG1c", "claw_tpGRN", "LB3a", "LB3b", "LB3c", "LB3d"]
CONTROL = ["JO-FV"]

neurons, edges = load_graph(DATA)
meta = neurons.set_index("bodyId")
jo = [t for t in find_types(neurons, [r"^JO-A", r"^JO-B"]) if not t.endswith("unclear")]
taste = gustatory_afferents(neurons) or [t for t in TASTE_FALLBACK if (neurons.type == t).any()]
print(f"gustatory afferent types ({'subclass' if gustatory_afferents(neurons) else 'FALLBACK list'}): {len(taste)}")
control = [t for t in CONTROL if (neurons.type == t).any()]
stim_all = list(LOOM_TUNING) + jo + taste + control
print(f"taste drive {taste}; control {control}; {len(jo)} JO types")


def gf_metrics(m, onsets, window=160):
    sp = m.spike_frame()
    gf = set(m.lif_ids[m.readout_index("DNp01")]); tt = set(m.lif_ids[m.readout_index("TTMn")])
    rows = []
    for on in onsets:
        w = sp[(sp.t_ms >= on) & (sp.t_ms < on + window)]
        g = w[w.bodyId.isin(gf)]; t = w[w.bodyId.isin(tt)]
        rows.append(dict(hit=int(len(g) > 0), per=(len(g) / g.bodyId.nunique()) if len(g) else np.nan,
                         gf=len(g), tt=len(t), lat=(g.t_ms.min() - on) if len(g) else np.nan))
    d = pd.DataFrame(rows)
    return dict(hit=d.hit.mean(), per=d.per.mean(), ratio=(d.tt.sum() / d.gf.sum()) if d.gf.sum() else np.nan,
                lat=d.lat.mean(), pop=m.population_rate_hz())


from flycns.graph import type_region
def mn9_metrics(m, onsets, verbose=False):
    r = readout_rates(m, onsets, ["MN9"], 250, meta)
    sp = m.spike_frame()
    reg = sp["type"].map(type_region)
    sc = sp["superclass"].fillna("")
    in_sez = (reg == "SEZ") | sc.isin(["cb_motor", "cb_sensory", "cb_sensory_tbc"])
    dur = m.t_ms / 1000.0
    out = dict(mn9=r.spikes_per_cell.mean() if len(r) else 0.0, pop=m.population_rate_hz(),
               pop_nonSEZ=(~in_sez).sum() / len(m.lif_ids) / dur)
    if verbose:
        mot = sp[sc == "cb_motor"].groupby("type").size().sort_values(ascending=False)
        n_mot = m.meta.loc[m.lif_ids][m.meta.loc[m.lif_ids, "superclass"] == "cb_motor"].groupby("type").size()
        if len(mot):
            print("  cb_motor types firing (spikes / n cells):", ", ".join(f"{t}:{v}/{n_mot.get(t, '?')}" for t, v in mot.head(8).items()),
                  f"  [{len(mot)} of {len(n_mot)} motor types active]")
    return out


rows = []
t_all = time.time()
for gsez in SEZ_GRID:
    gs, gr, gl = 1.0, 1.0, gsez        # columns kept for table compatibility; 'local' = SEZ gain
    t0 = time.time()
    p = LIFParams(region_gains={"SEZ": gsez, "other": 1.0})
    m = CNSModel(neurons, edges, stim_all, p, electrical=True)
    if not rows:
        print(f"SEZ definition: {getattr(m, 'region_source', 'unknown')}")
    m.store()
    # escape
    on = [o for o, _ in loom_protocol(m, n_trials=N)]; E = gf_metrics(m, on); m.restore()
    # JO alone
    on = pulse_protocol(m, {t: 150 for t in jo}, n_trials=N, pulse_ms=60); A1 = gf_metrics(m, on); m.restore()
    # JO + loom (loom tuning at gain 1 plus JO 150 Hz)
    tuning = m.stim["type"].map(LOOM_TUNING).fillna(0).values; is_jo = m.stim["type"].isin(jo).values
    on = []
    for _ in range(N):
        m.run(300); on.append(m.t_ms); m.set_stim_rates(tuning * 5.0 * 0.6 + is_jo * 150); m.run(60)
        m.set_stim_rates(np.zeros(len(m.stim)))
    m.run(300); A2 = gf_metrics(m, on); m.restore()
    # loom alone at fixed gain 1 for the A2 comparison
    on = []
    for _ in range(N):
        m.run(300); on.append(m.t_ms); m.set_stim_rates(tuning * 5.0 * 0.6); m.run(60); m.set_stim_rates(np.zeros(len(m.stim)))
    m.run(300); L0 = gf_metrics(m, on); m.restore()
    # feeding and specificity
    on = pulse_protocol(m, {t: 50 for t in taste}, n_trials=N); F = mn9_metrics(m, on, verbose=True)
    A.raster(m, on, ["MN9"], window_ms=250, bin_ms=8, max_trials=1); A.region_bar(m, on); m.restore()
    on = pulse_protocol(m, {t: 50 for t in control}, n_trials=N); S = mn9_metrics(m, on); m.restore()

    checks = dict(
        E1=(E["hit"] >= 0.8) and (1.0 <= (E["per"] or 0) <= 2.0),
        E2=0.8 <= (E["ratio"] if E["ratio"] == E["ratio"] else 0) <= 1.2,
        E3=E["pop"] < 0.05,
        A1=A1["hit"] <= 0.3,
        A2=(A2["hit"] > L0["hit"]) or (A2["lat"] < (L0["lat"] if L0["lat"] == L0["lat"] else 1e9) - 0.5)
           or (A2["per"] or 0) > (L0["per"] or 0),
        F1=F["mn9"] > 0.5,
        F2=F["pop"] < 0.05,
        S1=S["mn9"] < 0.2,
    )
    score = sum(checks.values())
    dev = abs(gsez - 1)
    row = dict(sensory=gs, relay=gr, local=gl, score=score, dev=dev, **{f"chk_{k}": v for k, v in checks.items()},
               gf_hit=E["hit"], gf_per=E["per"], ttmn_ratio=E["ratio"], pop_loom=E["pop"],
               jo_hit=A1["hit"], jo_loom_hit=A2["hit"], loom_hit=L0["hit"],
               mn9_taste=F["mn9"], pop_taste=F["pop"], pop_taste_nonSEZ=F["pop_nonSEZ"],
               mn9_control=S["mn9"], wall_s=round(time.time() - t0))
    rows.append(row)
    flags = "".join(k if v else "." for k, v in checks.items()).replace("E", "E").replace("A", "A")
    print(f"SEZ gain {gsez:<4} score {score}/8  "
          f"[{' '.join(k for k, v in checks.items() if v):<24}]  GF hit {E['hit']:.2f} per {E['per'] or 0:.2f} "
          f"| JO hit {A1['hit']:.2f} | MN9 taste {F['mn9']:.2f} ctrl {S['mn9']:.2f} pop {F['pop']:.3f} "
          f"(outside SEZ {F['pop_nonSEZ']:.4f}) | {row['wall_s']}s")
    pd.DataFrame(rows).to_csv(OUT / "grid.csv", index=False)

grid = pd.DataFrame(rows).sort_values(["score", "dev"], ascending=[False, True])
best = grid.iloc[0]
print(f"\nbest: SEZ gain {best.local}  score {int(best.score)}/8   "
      f"({(time.time()-t_all)/60:.0f} min total)")
print(grid[["sensory", "relay", "local", "score"] + [c for c in grid.columns if c.startswith("chk_")]].head(8).to_string(index=False))
(OUT / "best.json").write_text(json.dumps(dict(region_gains=dict(SEZ=float(best.local), other=1.0), score=int(best.score),
                                               checks={c[4:]: bool(best[c]) for c in grid.columns if c.startswith("chk_")}),
                                          indent=2))
