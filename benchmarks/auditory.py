"""
Benchmark 3: Johnston's organ -> giant fiber (mechanosensory input to escape).

Why this circuit: JON -> GF is a documented MIXED synapse that is PRIMARILY
ELECTRICAL (shakB); shakB2 mutants lose the JON-evoked GF current (Pezier &
Blagburn 2013; Yorozu et al. 2009). In vivo, sound alone leaves GF
subthreshold and summates with visual drive (von Reyn 2014). That yields
ablation predictions with opposite signs to the GF->TTMn case:

  A1 JO-A/B drive alone (pulse-song rates): GF response probability low, <= 0.3
  A2 REPORTED, NOT SCORED: effect of JO drive on GF response to a NEAR-THRESHOLD loom.
     The working point is calibrated per run (sweep of loom gain; take the largest gain
     with GF hit rate <= 0.3 and mean depolarisation >= 2 mV), because the direction of
     the effect depends on where the loom sits relative to threshold: with a loom that
     already fires GF there is only room to suppress, with one that never depolarises GF
     only room to facilitate. Earlier runs reported both directions at fixed gains (0.55
     -> 0.20 with the ordinal tuning at gain 0.6; 0.00 -> 0.07 with the measured tuning
     at the same gain) before this calibration existed; neither is quotable. The in vivo
     direction is not established in our references (von Reyn 2014 tested visual-visual
     integration; Pezier & Blagburn 2013 measured the JON-evoked GF potential alone).
  A3 pathway ablation: MaleCNS EM annotates the mixed JON-GF contact as
     chemical (679 synapses). Removing BOTH the electrical model and those
     EM edges must abolish the JO-evoked GF depolarisation (measured from GF
     membrane voltage, Pezier & Blagburn 2013); with the electrical model it
     must exceed 0.5 mV. "chem_only_EM" (EM as annotated, no electrical) is
     reported for reference: it double-counts nothing but mislabels the contact.
  A4 rewired null, chemical only: no JO -> GF depolarisation (the electrical
     pairs are declared anatomy and bypass rewiring, so they are removed here)
  A5 CNS stays quiet (< 0.05 Hz/neuron)

JO-A and JO-B are the vibration-sensitive subgroups (Kamikouchi et al. 2009);
pulse song drives them at ~100-200 Hz burst rates.
Writes results/auditory/{report.json, provenance.md, *.csv}
"""
import json, os, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flycns import load_graph, rewire_null, CNSModel, LIFParams, loom_protocol
from flycns.graph import drop_mixed_chemical
from flycns.protocols import LOOM_TUNING
from flycns.bench import find_types, pulse_protocol
from flycns.bench import provenance
from flycns import ascii as A

OUT = Path("results/auditory"); OUT.mkdir(parents=True, exist_ok=True)
N_TRIALS = int(sys.argv[1]) if len(sys.argv) > 1 else 40
# A2 compares two response probabilities, so it needs far more trials than the pass/fail
# checks: 20 trials cannot separate 0.15 from 0.30 (3 vs 6 events). Set FLYCNS_A2_TRIALS
# to run it properly (>=120 per arm); with fewer, A2 reports "unresolved".
A2_TRIALS = int(os.environ.get("FLYCNS_A2_TRIALS", str(N_TRIALS)))
DATA = sys.argv[2] if len(sys.argv) > 2 else "data"
JO_HZ = 150
SUB_LOOM_GAIN = None    # calibrated per run: largest gain with hit rate <= 0.3 and GF depolarised

neurons, edges = load_graph(DATA)
A.sketch("auditory")
meta = neurons.set_index("bodyId")
jo_types = [t for t in find_types(neurons, [r"^JO-A", r"^JO-B"]) if not t.endswith("unclear")]
print("JO vibration types:", jo_types)
if not jo_types:
    raise SystemExit("no JO-A/JO-B types found")
report = {"benchmark": "auditory", "n_trials": N_TRIALS, "jo_types": jo_types, "jo_hz": JO_HZ,
          "arms": {}, "checks": {}}


def gf_stats(model, onsets, window_ms=160):
    sp = model.spike_frame()
    gf_ids = set(model.lif_ids[model.readout_index("DNp01")])
    rows = []
    for on in onsets:
        w = sp[(sp.t_ms >= on) & (sp.t_ms < on + window_ms) & sp.bodyId.isin(gf_ids)]
        rows.append(dict(gf=len(w) / max(len(gf_ids), 1), hit=int(len(w) > 0),
                         lat=(w.t_ms.min() - on) if len(w) else np.nan))
    df = pd.DataFrame(rows)
    return dict(gf_per_cell=float(df.gf.mean()), gf_hit=float(df.hit.mean()),
                gf_lat_ms=float(df.lat.mean()), pop_rate_hz=float(model.population_rate_hz()))


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def two_prop_ci(k1, n1, k2, n2, z=1.96):
    """95% CI on p2 - p1 (Newcombe hybrid score interval)."""
    l1, u1 = wilson(k1, n1); l2, u2 = wilson(k2, n2)
    d = k2 / n2 - k1 / n1
    lo = d - z * np.sqrt((k1 / n1 - l1) ** 2 / z ** 2 * z ** 2 + (u2 - k2 / n2) ** 2 / z ** 2 * z ** 2)
    hi = d + z * np.sqrt((u1 - k1 / n1) ** 2 / z ** 2 * z ** 2 + (k2 / n2 - l2) ** 2 / z ** 2 * z ** 2)
    return d, lo, hi


def run(label, stim_types, electrical=True, edge_df=edges, jo=True, loom=False, n_trials=None):
    t0 = time.time()
    m = CNSModel(neurons, edge_df, stim_types, LIFParams(), electrical=electrical,
                 monitor_types=("DNp01",))
    N = N_TRIALS if n_trials is None else n_trials
    if loom and jo:
        # loom protocol drives LC types; add JO at fixed rate during each loom burst
        tuning = m.stim["type"].map(LOOM_TUNING).fillna(0).values
        is_jo = m.stim["type"].isin(jo_types).values
        onsets = []
        for _ in range(N):
            m.run(300); onsets.append(m.t_ms)
            m.set_stim_rates(tuning * 5.0 * SUB_LOOM_GAIN + is_jo * JO_HZ); m.run(60)
            m.set_stim_rates(np.zeros(len(m.stim)))
        m.run(300)
    elif loom:
        tuning = m.stim["type"].map(LOOM_TUNING).fillna(0).values
        onsets = []
        for _ in range(N):
            m.run(300); onsets.append(m.t_ms)
            m.set_stim_rates(tuning * 5.0 * SUB_LOOM_GAIN); m.run(60)
            m.set_stim_rates(np.zeros(len(m.stim)))
        m.run(300)
    else:
        onsets = pulse_protocol(m, {t: JO_HZ for t in stim_types}, n_trials=N, pulse_ms=60)
    s = gf_stats(m, onsets); s["wall_s"] = round(time.time() - t0, 1); s["n_trials"] = len(onsets)
    A.raster(m, onsets, ["DNp01"], window_ms=100, max_trials=2)
    s["gf_depol_mV"] = float(np.nanmean(m.peak_depolarization_mV(onsets, 100)))
    print(f"[{label}] GF {s['gf_per_cell']:.2f}/cell hit {s['gf_hit']:.2f} lat {s['gf_lat_ms']:.1f} ms "
          f"depol {s['gf_depol_mV']:.2f} mV | CNS {s['pop_rate_hz']:.4f} Hz | {s['wall_s']}s")
    report["arms"][label] = s
    return s


# ---- calibrate the near-threshold loom working point
CAL_GAINS = [0.4, 0.6, 0.8, 1.0, 1.3]
cal = []
for g in CAL_GAINS:
    SUB_LOOM_GAIN = g
    globals()["SUB_LOOM_GAIN"] = g
    c = run(f"cal_loom_gain_{g}", list(LOOM_TUNING), electrical=True, jo=False, loom=True)
    cal.append((g, c["gf_hit"], c["gf_depol_mV"]))
    print(f"  calibration: loom gain {g} -> GF hit {c['gf_hit']:.2f}, depol {c['gf_depol_mV']:.2f} mV")
ok = [(g, h, d) for g, h, d in cal if h <= 0.3 and d >= 2.0]
SUB_LOOM_GAIN = max(g for g, _, _ in ok) if ok else CAL_GAINS[0]
globals()["SUB_LOOM_GAIN"] = SUB_LOOM_GAIN
report["arms"]["A2_calibration"] = dict(gains=[dict(gain=g, hit=h, depol_mV=d) for g, h, d in cal],
                                        chosen_gain=SUB_LOOM_GAIN,
                                        rule="largest gain with GF hit <= 0.3 and depol >= 2 mV")
print(f"A2 working point: loom gain {SUB_LOOM_GAIN}")

a1 = run("A1_jo_alone_electrical", jo_types, electrical=True)
a3 = run("A3_jo_alone_no_electrical_no_mixed_edges", jo_types, electrical=False,
         edge_df=drop_mixed_chemical(edges, neurons))
ref = run("ref_jo_alone_chem_only_EM_as_annotated", jo_types, electrical=False)
lo = run("subthreshold_loom_alone", list(LOOM_TUNING), electrical=True, jo=False, loom=True, n_trials=A2_TRIALS)
a2 = run("A2_loom_plus_jo", list(LOOM_TUNING) + jo_types, electrical=True, jo=True, loom=True, n_trials=A2_TRIALS)
a4 = run("A4_jo_rewired_null_chem", jo_types, electrical=False,
         edge_df=rewire_null(drop_mixed_chemical(edges, neurons)))

report["checks"] = {
    "A1_jo_alone_mostly_subthreshold": a1["gf_hit"] <= 0.3,
    "A2_jo_effect_on_near_threshold_loom_REPORTED": None,
    "A3_declared_pairs_carry_jo_drive": (a1["gf_depol_mV"] > 0.5) and (a3["gf_depol_mV"] < 0.5 * a1["gf_depol_mV"]),
    "A4_null_silent": a4["gf_hit"] < 0.1 and a4["gf_depol_mV"] < 0.5,
    "A5_cns_quiet": max(a1["pop_rate_hz"], a2["pop_rate_hz"]) < 0.05,
}
k1, n1 = int(round(lo["gf_hit"] * lo["n_trials"])), lo["n_trials"]
k2, n2 = int(round(a2["gf_hit"] * a2["n_trials"])), a2["n_trials"]
d, dlo, dhi = two_prop_ci(k1, n1, k2, n2)
resolved = (dlo > 0) or (dhi < 0)
direction = ("facilitation" if d > 0 else "suppression" if d < 0 else "none") if resolved else "UNRESOLVED"
report["arms"]["A2_prediction"] = dict(loom_gain=SUB_LOOM_GAIN, n_trials_per_arm=n1,
    loom_alone_hit=lo["gf_hit"], loom_alone_ci95=list(wilson(k1, n1)),
    loom_plus_jo_hit=a2["gf_hit"], loom_plus_jo_ci95=list(wilson(k2, n2)),
    difference=d, difference_ci95=[dlo, dhi], resolved=bool(resolved), direction=direction,
    loom_alone_depol_mV=lo["gf_depol_mV"], loom_plus_jo_depol_mV=a2["gf_depol_mV"],
    note=("model measurement at a calibrated near-threshold working point; in vivo direction not "
          "established in cited references. With few trials the difference is not resolvable: "
          "set FLYCNS_A2_TRIALS>=120 for a usable estimate."))
print(f"A2 (reported, {n1} trials/arm): loom alone {lo['gf_hit']:.2f} -> with JO {a2['gf_hit']:.2f}; "
      f"difference {d:+.2f} (95% CI {dlo:+.2f} to {dhi:+.2f}) -> {direction}")
report["pass"] = all(v for v in report["checks"].values() if v is not None)
report["provenance"] = provenance()
(OUT / "report.json").write_text(json.dumps(report, indent=2, default=float))
prov = ["# Provenance: auditory (JON -> GF)", "", f"JO types: {jo_types} at {JO_HZ} Hz",
        f"A2 reported (not scored): JO drive {report['arms']['A2_prediction']['direction']} of near-threshold loom response",
        "Electrical JON->GF calibrated to a 3 mV compound GF potential at 150 Hz (Pezier & Blagburn 2013)",
        "", "## Checks"] + [f"- {k}: {'REPORTED' if v is None else ('PASS' if v else 'FAIL')}" for k, v in report["checks"].items()]
(OUT / "provenance.md").write_text("\n".join(prov))
print("\nCHECKS"); [print(f"  {'REPORTED' if v is None else ('PASS' if v else 'FAIL')}  {k}") for k, v in report["checks"].items()]
print(f"\nBENCHMARK {'PASS' if report['pass'] else 'FAIL'}  -> {OUT}/report.json")
