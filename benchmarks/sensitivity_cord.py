"""
Do the results depend on applying central-brain single-neuron parameters to the nerve cord?

Shiu et al. 2024's LIF values (tau_m 20 ms, threshold 7 mV above rest, refractory 2.2 ms) were
chosen for central-brain neurons and we apply them to every neuron, cord included. Motor neurons
are large, low-resistance cells; leg motor neuron recordings (Azevedo et al. 2020) put them in a
different regime. This sweeps cord parameters ONLY (superclasses vnc_motor, vnc_intrinsic,
vnc_sensory, vnc_efferent) and reports which checks move.

Reference points for the sweep:
  tau_m       5, 10, 20, 40 ms          (20 = the value in use)
  threshold   4, 7, 10 mV above rest    (7 = the value in use)
  refractory  1.0, 2.2, 5.0 ms          (2.2 = the value in use)
One parameter is varied at a time; everything else stays at the declared set.

Reads as: if the escape checks do not move, the result does not rest on cord parameters, which
is expected because the GF->TTMn relay is electrical and its spikelet dominates. If they move,
the affected checks need a caveat until cord parameters are fitted.

Run:  python benchmarks/sensitivity_cord.py [n_trials] [data_dir]
Writes results/sensitivity_cord/grid.csv
"""
import json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns import load_graph, CNSModel, LIFParams, loom_protocol
from flycns.protocols import LOOM_TUNING

N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
DATA = sys.argv[2] if len(sys.argv) > 2 else "data"
OUT = Path("results/sensitivity_cord"); OUT.mkdir(parents=True, exist_ok=True)
CORD = ["vnc_motor", "vnc_intrinsic", "vnc_sensory", "vnc_efferent"]
BASE = LIFParams()
SWEEP = ([("tau_m_ms", v) for v in (5.0, 10.0, 20.0, 40.0)]
         + [("v_thresh_mV", BASE.v_rest_mV + g) for g in (4.0, 7.0, 10.0)]
         + [("refractory_ms", v) for v in (1.0, 2.2, 5.0)])

neurons, edges = load_graph(DATA)
rows = []
for param, val in SWEEP:
    t0 = time.time()
    over = {sc: {param: val} for sc in CORD}
    m = CNSModel(neurons, edges, list(LOOM_TUNING), LIFParams(superclass_params=over), electrical=True)
    trials = loom_protocol(m, n_trials=N)
    sp = m.spike_frame()
    gf = set(m.lif_ids[m.readout_index("DNp01")]); tt = set(m.lif_ids[m.readout_index("TTMn")])
    per, hit, lat, ngf, ntt = [], [], [], 0, 0
    for on, _ in trials:
        w = sp[(sp.t_ms >= on) & (sp.t_ms < on + 160)]
        g = w[w.bodyId.isin(gf)]; t = w[w.bodyId.isin(tt)]
        ngf += len(g); ntt += len(t); hit.append(int(len(g) > 0))
        per.append(len(g) / g.bodyId.nunique() if len(g) else np.nan)
        if len(g) and len(t):
            lat.append(t.t_ms.min() - g.t_ms.min())
    n_cord = sum(x["n"] for x in m.superclass_params_applied)
    row = dict(param=param, value=val, n_cord_neurons=n_cord,
               gf_per_response=float(np.nanmean(per)) if np.isfinite(per).any() else 0.0,
               gf_hit=float(np.mean(hit)), ttmn_ratio=(ntt / ngf) if ngf else np.nan,
               ttmn_lag_ms=float(np.mean(lat)) if lat else np.nan,
               pop_rate_hz=float(m.population_rate_hz()), wall_s=round(time.time() - t0))
    row["E_spikes"] = 1.0 <= row["gf_per_response"] <= 2.0
    row["E_prob"] = 0.5 <= row["gf_hit"] <= 1.0
    row["E_relay"] = 0.8 <= (row["ttmn_ratio"] if np.isfinite(row["ttmn_ratio"]) else 0) <= 1.2
    row["E_lag"] = 0.5 <= (row["ttmn_lag_ms"] if np.isfinite(row["ttmn_lag_ms"]) else 0) <= 1.5
    row["E_quiet"] = row["pop_rate_hz"] < 0.05
    rows.append(row)
    flags = " ".join(k[2:] for k in ("E_spikes", "E_prob", "E_relay", "E_lag", "E_quiet") if row[k])
    print(f"cord {param:<14} = {val:>6}  ({n_cord:,} neurons) | GF {row['gf_per_response']:.2f}/resp "
          f"hit {row['gf_hit']:.2f} | TTMn/GF {row['ttmn_ratio']:.2f} lag {row['ttmn_lag_ms']:.2f} ms "
          f"| pop {row['pop_rate_hz']:.4f} | [{flags}] {row['wall_s']}s")
    pd.DataFrame(rows).to_csv(OUT / "grid.csv", index=False)

df = pd.DataFrame(rows)
checks = ["E_spikes", "E_prob", "E_relay", "E_lag", "E_quiet"]
moved = {c: int((~df[c]).sum()) for c in checks}
spread = {c: round(float(df[c.replace("E_spikes", "gf_per_response").replace("E_prob", "gf_hit")
                          .replace("E_relay", "ttmn_ratio").replace("E_lag", "ttmn_lag_ms")
                          .replace("E_quiet", "pop_rate_hz")].std()), 4) for c in checks}
summary = dict(n_settings=len(df), checks_failing_in_any_setting=moved, metric_std_across_sweep=spread,
               cord_superclasses=CORD,
               interpretation=("checks that never fail across a 8x range of cord tau_m, a 4-10 mV threshold "
                               "range and a 5x refractory range do not rest on cord single-neuron parameters"))
(OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
print("\nchecks failing in at least one setting:", moved)
print("metric spread across the sweep:", spread)
