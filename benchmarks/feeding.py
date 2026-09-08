"""
Benchmark 2: feeding (sugar GRN -> second-order SEZ -> MN9), a port of the
tests in Shiu et al. 2024 (Nature) from FlyWire to MaleCNS.

Criteria (Shiu 2024; Gordon & Scott 2009; McKellar 2020):
  F1 unilateral sugar GRN drive: contralateral MN9 > ipsilateral MN9
  F2 MN9 response is monotonic in sugar GRN rate over 10-200 Hz
  F3 bitter GRN co-activation suppresses MN9 relative to sugar alone
  F4 each of Fdg/Bract/Roundup/Zorro alone at 50 Hz is sufficient to drive MN9
  F5 rewired null: MN9 silent
Type discovery: sugar/bitter GRN names in MaleCNS are unverified; the script
lists candidates and falls back to second-order neurons (Shiu: Fdg sufficient).
Writes results/feeding/{report.json, provenance.md, *.csv}
"""
import json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns import load_graph, rewire_null, CNSModel, LIFParams
from flycns.graph import (FEEDING_SECOND_ORDER, FEEDING_MOTOR, SUGAR_GRN_PATTERNS,
                          BITTER_GRN_PATTERNS, TASTE_SENSORY_FALLBACK)
from flycns.bench import find_types, present_types, pulse_protocol, readout_rates
from flycns.graph import infer_side

OUT = Path("results/feeding"); OUT.mkdir(parents=True, exist_ok=True)
N_TRIALS = int(sys.argv[1]) if len(sys.argv) > 1 else 6
DATA = sys.argv[2] if len(sys.argv) > 2 else "data"
SUGAR = sys.argv[3].split(",") if len(sys.argv) > 3 else None   # override type names
BITTER = sys.argv[4].split(",") if len(sys.argv) > 4 else None

neurons, edges = load_graph(DATA)
meta = neurons.set_index("bodyId")

# ---------------------------------------------------------------- discovery
sugar_hits = find_types(neurons, SUGAR_GRN_PATTERNS)
bitter_hits = find_types(neurons, BITTER_GRN_PATTERNS)
so_have, so_miss = present_types(neurons, FEEDING_SECOND_ORDER)
# MaleCNS does not use Shiu's FlyWire names; discover MN9's strongest input types instead
_mn9 = meta.index[meta["type"] == "MN9"]
_exc = edges[edges.post.isin(_mn9) & (edges.sign > 0)].assign(t=lambda d: d.pre.map(meta["type"]))
mn9_inputs = _exc.groupby("t").weight.sum().sort_values(ascending=False)
if not so_have:
    so_have = list(mn9_inputs.head(4).index)
    print("second-order (Shiu names) absent; using MN9's top EXCITATORY input types:", so_have)
# sensory afferents feeding those excitatory inputs (superclass contains 'sensory')
_up = meta.index[meta["type"].isin(mn9_inputs.head(6).index)]
_s2 = edges[edges.post.isin(_up)].assign(t=lambda d: d.pre.map(meta["type"]),
                                        sc=lambda d: d.pre.map(meta["superclass"]).fillna(""))
taste_afferents = (_s2[_s2.sc.str.contains("sensory")].groupby("t").weight.sum()
                   .sort_values(ascending=False))
taste_afferents = list(taste_afferents[taste_afferents >= 50].index)
print("sensory afferents onto MN9 excitatory inputs (>=50 syn):", taste_afferents)
mn_have, mn_miss = present_types(neurons, FEEDING_MOTOR)
print("sugar GRN candidates:", sugar_hits or "NONE FOUND")
print("bitter GRN candidates:", bitter_hits or "NONE FOUND")
print("second-order present:", so_have, " missing:", so_miss)
print("motor present:", mn_have, " missing:", mn_miss)
if "MN9" not in mn_have:
    raise SystemExit("MN9 not found; cannot run feeding benchmark")

# drive priority: explicit argv list > sugar-specific GRN types if annotated >
# discovered MxLbN afferents onto MN9's excitatory inputs (GNG642, BM_Taste on v1.0)
specific = [t for t in sugar_hits if not t.startswith("BM_")]
sugar = SUGAR or specific or taste_afferents
bitter = BITTER or list(bitter_hits)
if not specific and not SUGAR:
    print(f"WARNING: no sugar-specific GRN type in this annotation; driving taste-region afferents "
          f"{sugar} (F3 bitter test skipped; F1/F2/F5 test taste -> MN9 wiring)")
if not sugar:
    raise SystemExit("no usable sensory drive; supply type names as argv[3]")

report = {"benchmark": "feeding", "n_trials": N_TRIALS, "sugar_types": sugar, "bitter_types": bitter,
          "second_order_present": so_have, "motor_present": mn_have, "arms": {}, "checks": {}}
WINDOW = 250


def mn9(model, onsets):
    df = readout_rates(model, onsets, ["MN9"], WINDOW, meta)
    return df


def build(stim_types, edge_df=edges):
    return CNSModel(neurons, edge_df, stim_types, LIFParams(), electrical=True)


# ---------------------------------------------------------------- F1 laterality
t0 = time.time()
sugar_ids = neurons[neurons.type.isin(sugar)].copy()
sugar_ids["side"] = infer_side(meta).loc[sugar_ids.bodyId].values
left_sugar = sugar_ids[sugar_ids.side == "L"]
print(f"sugar drive: {len(sugar_ids)} neurons, {len(left_sugar)} inferred left")
m = build(sugar)
rates = pd.Series(0.0, index=m.stim.bodyId.values)
rates[left_sugar.bodyId.values] = 50.0
onsets = []
for _ in range(N_TRIALS):
    m.run(300); onsets.append(m.t_ms); m.set_stim_rates(rates.loc[m.stim.bodyId].values); m.run(200)
    m.set_stim_rates(np.zeros(len(m.stim)))
m.run(300)
lat = mn9(m, onsets)
ipsi = lat[lat.side == "L"].spikes_per_cell.mean(); contra = lat[lat.side == "R"].spikes_per_cell.mean()
report["arms"]["F1_left_sugar_50Hz"] = dict(MN9_ipsi_L=float(ipsi), MN9_contra_R=float(contra),
                                            pop_rate_hz=float(m.population_rate_hz()), wall_s=round(time.time()-t0, 1))
report["checks"]["F1_contra_gt_ipsi"] = bool(contra > ipsi > -1 and contra > 0)
print(f"F1 left sugar 50 Hz: MN9 ipsi {ipsi:.2f}, contra {contra:.2f} spikes/cell  [{time.time()-t0:.0f}s]")
lat.to_csv(OUT / "F1_laterality.csv", index=False)

# ---------------------------------------------------------------- F2 dose response
curve = []
for hz in [10, 25, 50, 100, 200]:
    m = build(sugar)
    on = pulse_protocol(m, {t: hz for t in sugar}, n_trials=N_TRIALS)
    v = mn9(m, on).spikes_per_cell.mean()
    curve.append(dict(sugar_hz=hz, MN9_spikes_per_cell=float(v), pop_rate_hz=float(m.population_rate_hz())))
    print(f"F2 sugar {hz:>3} Hz -> MN9 {v:.2f} spikes/cell, CNS {m.population_rate_hz():.4f} Hz")
curve = pd.DataFrame(curve); curve.to_csv(OUT / "F2_dose.csv", index=False)
report["arms"]["F2_dose"] = curve.to_dict("records")
diffs = np.diff(curve.MN9_spikes_per_cell.values)
report["checks"]["F2_monotonic"] = bool((diffs >= -1e-9).all() and curve.MN9_spikes_per_cell.iloc[-1] > 0)

# ---------------------------------------------------------------- F3 bitter suppression
if bitter:
    m = build(sugar + bitter)
    on = pulse_protocol(m, {**{t: 50 for t in sugar}, **{t: 50 for t in bitter}}, n_trials=N_TRIALS)
    both = mn9(m, on).spikes_per_cell.mean()
    alone = float(curve.loc[curve.sugar_hz == 50, "MN9_spikes_per_cell"].iloc[0])
    report["arms"]["F3_bitter"] = dict(MN9_sugar_alone=alone, MN9_sugar_plus_bitter=float(both))
    report["checks"]["F3_bitter_suppresses"] = bool(both < alone)
    print(f"F3 sugar alone {alone:.2f} vs sugar+bitter {both:.2f}")
else:
    report["checks"]["F3_bitter_suppresses"] = None
    print("F3 skipped: no bitter GRN type resolved")

# ---------------------------------------------------------------- F4 second-order sufficiency
suff = {}
for typ in [t for t in (["Fdg", "Bract", "Roundup", "Zorro"] if "Fdg" in so_have else so_have) if t in so_have]:
    m = build([typ])
    on = pulse_protocol(m, {typ: 50}, n_trials=N_TRIALS)
    suff[typ] = float(mn9(m, on).spikes_per_cell.mean())
    print(f"F4 {typ} 50 Hz -> MN9 {suff[typ]:.2f}")
report["arms"]["F4_sufficiency"] = suff
report["checks"]["F4_second_order_sufficient"] = (bool(suff) and all(v > 0 for v in suff.values())) if suff else None

# ---------------------------------------------------------------- F5 null
m = build(sugar, rewire_null(edges))
on = pulse_protocol(m, {t: 50 for t in sugar}, n_trials=N_TRIALS)
nullv = float(mn9(m, on).spikes_per_cell.mean())
report["arms"]["F5_null"] = dict(MN9=nullv)
report["checks"]["F5_null_silent"] = bool(nullv < 0.1)
print(f"F5 rewired null -> MN9 {nullv:.2f}")

report["pass"] = all(v for v in report["checks"].values() if v is not None)
(OUT / "report.json").write_text(json.dumps(report, indent=2, default=float))
prov = ["# Provenance: feeding (port of Shiu et al. 2024 tests to MaleCNS)", "",
        f"sugar types used: {sugar}", f"bitter types used: {bitter or 'none resolved'}",
        f"second-order present: {so_have}", f"motor present: {mn_have}", "", "## Checks"]
prov += [f"- {k}: {'SKIP' if v is None else ('PASS' if v else 'FAIL')}" for k, v in report["checks"].items()]
(OUT / "provenance.md").write_text("\n".join(prov))
print("\nCHECKS"); [print(f"  {'SKIP' if v is None else ('PASS' if v else 'FAIL')}  {k}") for k, v in report["checks"].items()]
print(f"\nBENCHMARK {'PASS' if report['pass'] else 'FAIL'}  -> {OUT}/report.json")
