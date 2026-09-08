"""
Benchmark 2: feeding (sugar GRN -> second-order SEZ -> MN9), a port of the
tests in Shiu et al. 2024 (Nature) from FlyWire to MaleCNS.

Two declared drive arms, because MaleCNS does not annotate taste modality:
  ALL     every anatomically gustatory afferent type (labellar bristle, taste peg,
          pharyngeal sensillum; 42 types). Sensilla house sugar, bitter, water, salt
          and mechanosensory neurons together, so this arm is "every tastant at once"
          and is NOT expected to drive PER. Reported, not scored.
  SUBSET  the screen-selected MN9-driving gustatory types (GUSTATORY_MN9_DRIVING),
          the closest sugar-GRN proxy available; scored on F1/F2. Wiring-selected,
          not a modality annotation.
Criteria (Shiu 2024; Gordon & Scott 2009; McKellar 2020):
  F1 unilateral drive: MN9 laterality has the same sign as the structural 2-hop
     prediction (dynamics reproduce wiring), at the lowest dose-curve rate giving MN9 >= 2
     spikes/cell bilaterally. Shiu 2024's contralateral bias was for
     labellar SUGAR GRNs, which MaleCNS does not annotate: recorded as N/A, not scored.
  F2 MN9 response is monotonic in sugar GRN rate over 10-200 Hz
  F3 bitter GRN co-activation suppresses MN9 relative to sugar alone
  F4 each of Fdg/Bract/Roundup/Zorro alone at 50 Hz is sufficient to drive MN9
  F5 rewired null: MN9 silent
  F2b activity outside the SEZ stays < 0.02 Hz/neuron during taste drive
  F2c recruited motor neurons fire < 100 Hz on average (Azevedo 2020; McKellar 2020)
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
from flycns.graph import infer_side, gustatory_afferents, GUSTATORY_MN9_DRIVING
from flycns import ascii as A

OUT = Path("results/feeding"); OUT.mkdir(parents=True, exist_ok=True)
N_TRIALS = int(sys.argv[1]) if len(sys.argv) > 1 else 6
DATA = sys.argv[2] if len(sys.argv) > 2 else "data"
SUGAR = sys.argv[3].split(",") if len(sys.argv) > 3 else None   # override type names
BITTER = sys.argv[4].split(",") if len(sys.argv) > 4 else None

neurons, edges = load_graph(DATA)
A.sketch("feeding")
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
gust = gustatory_afferents(neurons)
subset = [t for t in GUSTATORY_MN9_DRIVING if (neurons.type == t).any()]
sugar = SUGAR or specific or subset or gust or taste_afferents
report_all_arm = gust if (gust and not SUGAR and not specific) else None
if gust and not SUGAR and not specific:
    print(f"gustatory afferents by subclass (labellar bristle, taste peg, pharyngeal sensillum): {len(gust)} types")
bitter = BITTER or list(bitter_hits)
if not specific and not SUGAR:
    print(f"WARNING: no sugar-specific GRN type in this annotation; driving taste-region afferents "
          f"{sugar} (F3 bitter test skipped; F1/F2/F5 test taste -> MN9 wiring)")
if not sugar:
    raise SystemExit("no usable sensory drive; supply type names as argv[3]")

report = {"benchmark": "feeding", "n_trials": N_TRIALS, "sugar_types": sugar, "bitter_types": bitter,
          "drive_arm": "explicit" if SUGAR else "annotated_sugar" if specific else "screen_selected_subset" if subset else "all_gustatory",
          "second_order_present": so_have, "motor_present": mn_have, "arms": {}, "checks": {}}
WINDOW = 250


def mn9(model, onsets):
    df = readout_rates(model, onsets, ["MN9"], WINDOW, meta)
    return df


def build(stim_types, edge_df=edges):
    return CNSModel(neurons, edge_df, stim_types, LIFParams(), electrical=True)


# ---------------------------------------------------------------- F2 dose response
curve = []
for hz in [10, 25, 50, 100, 200]:
    m = build(sugar)
    on = pulse_protocol(m, {t: hz for t in sugar}, n_trials=N_TRIALS)
    v = mn9(m, on).spikes_per_cell.mean()
    curve.append(dict(sugar_hz=hz, MN9_spikes_per_cell=float(v), pop_rate_hz=float(m.population_rate_hz())))
    if hz == 50:
        A.raster(m, on, ["MN9"], window_ms=250, bin_ms=8, max_trials=2); A.region_bar(m, on)
    print(f"F2 sugar {hz:>3} Hz -> MN9 {v:.2f} spikes/cell, CNS {m.population_rate_hz():.4f} Hz")
    if hz == 50:
        from flycns.graph import type_region
        spf = m.spike_frame(); scc = spf["superclass"].fillna("")
        in_sez = (spf["type"].map(type_region) == "SEZ") | scc.isin(["cb_motor", "cb_sensory", "cb_sensory_tbc"])
        dur = m.t_ms / 1000.0
        report["arms"]["F2b_nonSEZ_rate_hz"] = float((~in_sez).sum() / len(m.lif_ids) / dur)
        mot = spf[scc == "cb_motor"]
        stim_s = N_TRIALS * 0.25
        per_cell = mot.groupby("bodyId").size() / stim_s
        report["arms"]["F2c_motor_mean_hz"] = float(per_cell.mean()) if len(per_cell) else 0.0
        n_types = m.meta.loc[m.lif_ids][m.meta.loc[m.lif_ids, "superclass"] == "cb_motor"].type.nunique()
        report["arms"]["motor_types_active"] = {t: int(v) for t, v in mot.groupby("type").size().sort_values(ascending=False).head(12).items()}
        print(f"  outside-SEZ rate {report['arms']['F2b_nonSEZ_rate_hz']:.4f} Hz | recruited motor mean {report['arms']['F2c_motor_mean_hz']:.1f} Hz "
              f"| {mot.type.nunique()} of {n_types} cb_motor types active: {list(report['arms']['motor_types_active'])[:8]}")
        report["checks"]["F2b_extra_SEZ_quiet"] = bool(report["arms"]["F2b_nonSEZ_rate_hz"] < 0.02)
        report["checks"]["F2c_motor_rate_physiological"] = bool(0 < report["arms"]["F2c_motor_mean_hz"] < 100)
curve = pd.DataFrame(curve); curve.to_csv(OUT / "F2_dose.csv", index=False)
if report_all_arm:
    m = build(report_all_arm)
    on = pulse_protocol(m, {t: 50 for t in report_all_arm}, n_trials=N_TRIALS)
    v_all = float(mn9(m, on).spikes_per_cell.mean())
    report["arms"]["ALL_gustatory_50Hz"] = dict(n_types=len(report_all_arm), MN9_spikes_per_cell=v_all,
                                              pop_rate_hz=float(m.population_rate_hz()))
    print(f"ALL gustatory arm ({len(report_all_arm)} types, every tastant at once) at 50 Hz -> MN9 {v_all:.2f} "
          f"spikes/cell (reported, not scored; expected ~0)")
report["arms"]["F2_dose"] = curve.to_dict("records")
diffs = np.diff(curve.MN9_spikes_per_cell.values)
report["checks"]["F2_monotonic"] = bool((diffs >= -1e-9).all() and curve.MN9_spikes_per_cell.iloc[-1] > 0)

# ---------------------------------------------------------------- F1 laterality
t0 = time.time()
sugar_ids = neurons[neurons.type.isin(sugar)].copy()
sugar_ids["side"] = infer_side(meta).loc[sugar_ids.bodyId].values
left_sugar = sugar_ids[sugar_ids.side == "L"]
print(f"sugar drive: {len(sugar_ids)} neurons, {len(left_sugar)} inferred left")
m = build(sugar)
rates = pd.Series(0.0, index=m.stim.bodyId.values)
# rate: lowest dose-curve rate at which bilateral drive gave MN9 >= 1 spike/cell,
# so a left/right ratio is measurable (unilateral drive is half the input)
_ok = curve[curve.MN9_spikes_per_cell >= 2.0]   # >= 2 spikes/cell so a side ratio is measurable
F1_HZ = float(_ok.sugar_hz.iloc[0]) if len(_ok) else 100.0
rates[left_sugar.bodyId.values] = F1_HZ
onsets = []
for _ in range(N_TRIALS):
    m.run(300); onsets.append(m.t_ms); m.set_stim_rates(rates.loc[m.stim.bodyId].values); m.run(200)
    m.set_stim_rates(np.zeros(len(m.stim)))
m.run(300)
lat = mn9(m, onsets)
ipsi = lat[lat.side == "L"].spikes_per_cell.mean(); contra = lat[lat.side == "R"].spikes_per_cell.mean()
report["arms"]["F1_left_sugar"] = dict(rate_hz=F1_HZ, MN9_ipsi_L=float(ipsi), MN9_contra_R=float(contra),
                                            pop_rate_hz=float(m.population_rate_hz()), wall_s=round(time.time()-t0, 1))
# structural prediction: 2-hop synapse-weighted drive from the driven (left) afferents onto each MN9
_L = list(left_sugar.bodyId); _s1 = edges[edges.pre.isin(_L)]; _mid = _s1.groupby("post").weight.sum()
_s2 = edges[edges.pre.isin(_mid.index) & edges.post.isin(_mn9)].assign(w=lambda d: d.weight * d.pre.map(_mid))
_side = infer_side(meta); _struct = _s2.groupby(_s2.post.map(_side)).w.sum()
struct_L, struct_R = float(_struct.get("L", 0)), float(_struct.get("R", 0))
report["arms"]["F1_left_sugar"]["structural_2hop_L"] = struct_L
report["arms"]["F1_left_sugar"]["structural_2hop_R"] = struct_R
func_sign = np.sign(contra - ipsi); struct_sign = np.sign(struct_R - struct_L)
report["checks"]["F1_laterality_matches_wiring"] = bool((ipsi + contra) > 0 and func_sign == struct_sign)
# Shiu 2024 predicted contralateral > ipsilateral MN9 for LABELLAR SUGAR GRNs. MaleCNS v1.0
# has no sugar annotation; the driven subset is pharyngeal/taste-peg/LB3 dominated, so the
# comparison is not applicable and is recorded as such rather than scored.
report["checks"]["F1_shiu_contralateral_bias"] = None
report["arms"]["F1_left_sugar"]["note"] = ("Shiu 2024 contralateral prediction concerns labellar sugar GRNs; "
                                          "not testable without a modality annotation")
print(f"   structural 2-hop drive from left afferents: MN9-L {struct_L:,.0f} vs MN9-R {struct_R:,.0f}; "
      f"functional {'ipsi' if ipsi > contra else 'contra'} dominant -> "
      f"{'matches' if report['checks']['F1_laterality_matches_wiring'] else 'does NOT match'} wiring")
print(f"F1 left sugar {F1_HZ:.0f} Hz: MN9 ipsi {ipsi:.2f}, contra {contra:.2f} spikes/cell  [{time.time()-t0:.0f}s]")
lat.to_csv(OUT / "F1_laterality.csv", index=False)

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
