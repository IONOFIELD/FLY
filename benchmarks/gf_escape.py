"""
Benchmark 1: giant fiber escape circuit (loom -> LC4/LPLC2 -> GF -> TTMn/PSI).

Arms:  real+electrical | real chemical-only | rewired null (+electrical)
Pass criteria (published physiology):
  GF 1-2 spikes per responding trial, hit rate >= 0.8, latency 10-60 ms
  TTMn spikes == GF spikes (one-to-one), lag 0.5-1.5 ms
  whole-CNS rate < 0.05 Hz/neuron ; null GF hit rate < 0.1
Writes: results/gf_escape/{arm}.csv, spikes_real.parquet, report.json, provenance.md
"""
import json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns import load_graph, rewire_null, CNSModel, LIFParams, loom_protocol
from flycns.protocols import LOOM_TUNING
from flycns.graph import ELECTRICAL_SYNAPSES, INTRINSIC_OVERRIDES, SIGN_MAP
from flycns import ascii as A

OUT = Path("results/gf_escape"); OUT.mkdir(parents=True, exist_ok=True)
N_TRIALS = int(sys.argv[1]) if len(sys.argv) > 1 else 20
DATA = sys.argv[2] if len(sys.argv) > 2 else "data"


def score(model, trials, burst_ms=60, window_ms=100):
    sp = model.spike_frame()
    ro = {r: set(model.lif_ids[model.readout_index(r)]) for r in ["DNp01", "TTMn", "PSI"]}
    rows = []
    for k, (on, gain) in enumerate(trials):
        w = sp[(sp.t_ms >= on) & (sp.t_ms < on + burst_ms + window_ms)]
        gf = w[w.bodyId.isin(ro["DNp01"])]; tt = w[w.bodyId.isin(ro["TTMn"])]; ps = w[w.bodyId.isin(ro["PSI"])]
        rows.append(dict(
            trial=k, gain=gain,
            gf=len(gf) / max(len(ro["DNp01"]), 1), gf_lat=gf.t_ms.min() - on if len(gf) else np.nan,
            gf_per_firing_cell=(len(gf) / gf.bodyId.nunique()) if len(gf) else np.nan,
            ttmn=len(tt) / max(len(ro["TTMn"]), 1),
            ttmn_lag=tt.t_ms.min() - gf.t_ms.min() if (len(gf) and len(tt)) else np.nan,
            psi=len(ps) / max(len(ro["PSI"]), 1)))
    return pd.DataFrame(rows), sp


def run_arm(name, neurons, edges, electrical):
    t0 = time.time()
    m = CNSModel(neurons, edges, list(LOOM_TUNING), LIFParams(), electrical=electrical)
    trials = loom_protocol(m, n_trials=N_TRIALS)
    df, sp = score(m, trials)
    A.raster(m, [t for t, _ in trials], ["DNp01", "TTMn"], meta=m.meta)
    A.region_bar(m, [t for t, _ in trials], window_ms=160)
    summ = dict(arm=name, n_lif=int(len(m.lif_ids)), n_chem_syn=int(m.n_chem),
                n_stim=int(len(m.stim)), pop_rate_hz=float(m.population_rate_hz()),
                gf_per_cell=float(df.gf.mean()), gf_hit=float((df.gf > 0).mean()),
                gf_per_hit=float(df.gf_per_firing_cell.mean()) if (df.gf > 0).any() else 0.0,
                gf_lat_ms=float(df.gf_lat.mean()), ttmn_per_cell=float(df.ttmn.mean()),
                ttmn_lag_ms=float(df.ttmn_lag.mean()), psi_per_cell=float(df.psi.mean()),
                ttmn_to_gf_ratio=float(df.ttmn.sum() / df.gf.sum()) if df.gf.sum() else float("nan"),
                wall_s=round(time.time() - t0, 1))
    print(f"[{name}] GF {summ['gf_per_cell']:.2f}/cell hit {summ['gf_hit']:.2f} lat {summ['gf_lat_ms']:.1f} ms | "
          f"TTMn {summ['ttmn_per_cell']:.2f}/cell lag {summ['ttmn_lag_ms']:.2f} ms | "
          f"pop {summ['pop_rate_hz']:.4f} Hz | {summ['wall_s']}s")
    df.to_csv(OUT / f"{name}.csv", index=False)
    return summ, sp, m


neurons, edges = load_graph(DATA)
A.sketch("gf_escape")
report = {"benchmark": "gf_escape", "n_trials": N_TRIALS, "params": LIFParams().__dict__, "arms": {}}

s, sp, m = run_arm("real_electrical", neurons, edges, True)
report["arms"]["real_electrical"] = s
sp.to_parquet(OUT / "spikes_real.parquet", index=False)
s2, _, _ = run_arm("real_chemical_only", neurons, edges, False)
report["arms"]["real_chemical_only"] = s2
s3, _, _ = run_arm("rewired_null", neurons, rewire_null(edges), True)
report["arms"]["rewired_null"] = s3

r, c, n = s, s2, s3
checks = {
    # spikes per GF that fired (von Reyn 2014 recorded single GFs): looms often
    # leave GF subthreshold, so failures are scored by gf_hit, count by responders
    "gf_spikes_per_response_1_to_2": 1.0 <= r["gf_per_hit"] <= 2.0,
    "gf_hit_ge_0.8": r["gf_hit"] >= 0.8,
    "gf_latency_10_60ms": 10 <= r["gf_lat_ms"] <= 60,
    "ttmn_one_to_one": 0.8 <= r["ttmn_to_gf_ratio"] <= 1.2,
    "ttmn_lag_0.5_1.5ms": 0.5 <= r["ttmn_lag_ms"] <= 1.5,
    "cns_quiet": r["pop_rate_hz"] < 0.05,
    "electrical_necessary": c["ttmn_per_cell"] < 0.1 * max(r["ttmn_per_cell"], 1e-9),
    "null_silent": n["gf_hit"] < 0.1,
}
report["checks"] = checks
report["pass"] = all(checks.values())
(OUT / "report.json").write_text(json.dumps(report, indent=2, default=float))

prov = ["# Provenance: gf_escape", "", "## Neurotransmitter sign map", str(SIGN_MAP), "",
        "## Electrical synapses added (invisible to EM)"]
prov += [f"- {a} -> {b} ({'ipsilateral' if i else 'any side'}, spikelet {k}): {s}" for a, b, i, k, s in ELECTRICAL_SYNAPSES]
prov += ["", "## Intrinsic overrides"] + [f"- {t}.{p} = {v}: {s}" for t, p, v, s in INTRINSIC_OVERRIDES]
prov += ["", f"## Regional gain: SEZ x{LIFParams().region_gains['SEZ']} (fitted, benchmarks/fit_regional.py), other x1"]
prov += ["", "## Loom tuning (PLACEHOLDER ordinal, Turner et al. 2022 Fig 3A)", str(LOOM_TUNING), "",
         "## Checks"] + [f"- {k}: {'PASS' if v else 'FAIL'}" for k, v in checks.items()]
(OUT / "provenance.md").write_text("\n".join(prov))

print("\nCHECKS"); [print(f"  {'PASS' if v else 'FAIL'}  {k}") for k, v in checks.items()]
print(f"\nBENCHMARK {'PASS' if report['pass'] else 'FAIL'}  -> {OUT}/report.json")
