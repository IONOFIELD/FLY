"""
Does a tonic baseline open the taste -> MN9 pathway (disinhibition), and at what cost?

Signed path products from gustatory afferents to MN9 in MaleCNS are negative at two hops
and positive at three to five: the pathway is disinhibitory. A network with no spontaneous
activity cannot express that (Shiu et al. 2024 ran at 0 Hz and noted the same limitation).
This sweeps the tonic background rate with the regional gain OFF, and scores:

  F  taste drive (screen-selected gustatory subset, 50 Hz) -> MN9 spikes/cell   want > 0.5
  S  wind control (JO-FV, 50 Hz) -> MN9                                          want < 0.2
  E  loom -> GF spikes per responding cell and response probability             want 1-2, 0.5-1
  Q  spontaneous whole-CNS rate with no stimulus                     want 0.5-5 Hz per neuron
     (in vivo fly neurons are not silent; this is the point of the baseline, and also its cost)

Run:  python benchmarks/fit_baseline.py [n_trials] [data_dir]
Writes results/fit_baseline/grid.csv
"""
import json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns import load_graph, CNSModel, LIFParams, loom_protocol
from flycns.protocols import LOOM_TUNING
from flycns.bench import pulse_protocol, readout_rates
from flycns.graph import GUSTATORY_MN9_DRIVING

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
DATA = sys.argv[2] if len(sys.argv) > 2 else "data"
OUT = Path("results/fit_baseline"); OUT.mkdir(parents=True, exist_ok=True)
# (tonic depolarisation mV, membrane noise mV). Rest sits 7 mV below threshold, so these
# set how close to threshold cells idle and how much they fluctuate. Target: 0.5-5 Hz
# spontaneous, the in vivo range.
GRID = [(0.0, 0.0), (3.0, 1.0), (4.0, 1.5), (5.0, 1.5), (5.0, 2.5), (6.0, 2.0)]

neurons, edges = load_graph(DATA)
meta = neurons.set_index("bodyId")
taste = [t for t in GUSTATORY_MN9_DRIVING if (neurons.type == t).any()]
control = ["JO-FV"] if (neurons.type == "JO-FV").any() else []
stim_all = list(LOOM_TUNING) + taste + control
rows = []
for depol, noise in GRID:
    b = f"{depol}/{noise}"
    t0 = time.time()
    p = LIFParams(baseline_depol_mV=depol, baseline_noise_mV=noise,
                  region_gains={"SEZ": 1.0, "other": 1.0})                  # regional gain OFF
    m = CNSModel(neurons, edges, stim_all, p, electrical=True)
    m.store()
    m.run(500); quiet = m.population_rate_hz(); m.restore()                 # cost of baseline alone
    on = pulse_protocol(m, {t: 50 for t in taste}, n_trials=N)
    r = readout_rates(m, on, ["MN9"], 250, meta); mn9 = r.spikes_per_cell.mean() if len(r) else 0.0
    pop_taste = m.population_rate_hz(); m.restore()
    ctrl = 0.0
    if control:
        on = pulse_protocol(m, {t: 50 for t in control}, n_trials=N)
        rc = readout_rates(m, on, ["MN9"], 250, meta); ctrl = rc.spikes_per_cell.mean() if len(rc) else 0.0
        m.restore()
    tr = loom_protocol(m, n_trials=N)
    sp = m.spike_frame(); gf = set(m.lif_ids[m.readout_index("DNp01")])
    per, hit = [], []
    for on_ms, _ in tr:
        w = sp[(sp.t_ms >= on_ms) & (sp.t_ms < on_ms + 160) & sp.bodyId.isin(gf)]
        hit.append(int(len(w) > 0)); per.append(len(w) / w.bodyId.nunique() if len(w) else np.nan)
    row = dict(depol_mV=depol, noise_mV=noise, quiet_rate_hz=quiet, mn9_taste=mn9, mn9_control=ctrl,
               pop_taste=pop_taste, gf_per_response=np.nanmean(per) if np.isfinite(per).any() else 0.0,
               gf_hit=float(np.mean(hit)), wall_s=round(time.time() - t0))
    row["F"] = row["mn9_taste"] > 0.5; row["S"] = row["mn9_control"] < 0.2
    row["E"] = (1.0 <= row["gf_per_response"] <= 2.0) and (0.5 <= row["gf_hit"] <= 1.0)
    row["Q"] = 0.5 <= row["quiet_rate_hz"] <= 5.0
    row["ok"] = bool(row["F"] and row["S"] and row["E"] and row["Q"])
    rows.append(row)
    print(f"depol {depol:>4} mV noise {noise:>4} mV | spont {quiet:6.3f} Hz | MN9 taste {mn9:6.2f} ctrl {ctrl:5.2f} | "
          f"GF {row['gf_per_response']:.2f}/resp hit {row['gf_hit']:.2f} | "
          f"[{' '.join(k for k in 'FSEQ' if row[k])}] {'OK' if row['ok'] else ''} {row['wall_s']}s")
    pd.DataFrame(rows).to_csv(OUT / "grid.csv", index=False)
df = pd.DataFrame(rows); ok = df[df.ok]
print("\nsatisfying F,S,E,Q:", ok[["depol_mV","noise_mV","quiet_rate_hz","mn9_taste"]].to_string(index=False) if len(ok) else "none")
(OUT / "best.json").write_text(json.dumps(dict(
    baseline_depol_mV=(float(ok.iloc[0].depol_mV) if len(ok) else None),
    baseline_noise_mV=(float(ok.iloc[0].noise_mV) if len(ok) else None),
    rule="smallest tonic rate with taste->MN9 > 0.5, wind control < 0.2, escape intact, spontaneous rate 0.5-5 Hz",
    rationale="taste->MN9 is disinhibitory (signed path products negative at 2 hops, positive at 3-5); "
              "a silent network cannot express disinhibition (cf. Shiu et al. 2024 limitation)"), indent=2))
